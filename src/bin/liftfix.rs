/*!

# What is it

Fix annotations (tbl) based on VCF generated after running frame_stats.
It uses as input the FASTA produced by bcftools consensus.
For instance, bcftools consensus -f examples/input_diamond_sample_genome.fa -c out.chain out.vcf.gz > out_consensus.fa
The algo is basically:
    1) It reads annotation features from tbl;
    2) liftover features coordinates;
    3) Join features that were detected as interrupted ORFs in frame_stats
    4) Outputs features in GFF. The GFF have specs according to NCBI's annnotation.


### Usage, option summary and outputs

```text
USAGE:
    liftfix [OPTIONS] <INPUT>

ARGS:
    <INPUT>    Input tbl

OPTIONS:
    -v, --variants <VCF>       VCF
    -g, --genes <GENES>        Fasta file with genes sequences (DNA). The gene id must match,
                               locus_tag in tbl input.
    -a, --help                 Print help information
    -r <GENOME>                FASTA produced by bcftools consensus. Requires faidx index.
    -V, --version              Print version information
```

## Examples

```text
$ /target/debug/frame_stats -d examples/output_diamond_sample.txt
    -g examples/input_diamond_sample_cds.fa -r examples/input_diamond_sample_genome.fa
     examples/input_diamond_sample.tbl
Number of CDS = 100

$
```

*/

extern crate libflate;
extern crate uniprot;
extern crate ureq;

use bio::data_structures::annot_map::AnnotMap;
use bio::io::fasta::{self, IndexedReader, Record};
use bio_types::annot::pos::Pos;
use bio_types::strand::ReqStrand;
use noodles::bgzf::{self as bgzf};
use regex::Regex;
use rs_genomesub::liftover::{LiftBlock, LiftBlockBuild, LiftOver};
use rs_genomesub::tbl::FeatureType;
use std::cell::RefCell;
use std::error::Error;
use std::io::prelude::*;
use std::io::BufReader;
use std::{io, u64};

use crate::rs_genomesub::chain;
use crate::rs_genomesub::liftover;

// use nom::error::dbg_dmp;
// use rs_genomesub::tbl::Feature;
// use uniprot::uniprot::DbReference;
use ureq::Response;

extern crate rs_genomesub;
use crate::rs_genomesub::blast_tabular::*;
use crate::rs_genomesub::tbl;

use itertools::Itertools;
use std::collections::HashMap;

use std::fs;
use std::fs::File;
// use std::io::Write;
use std::path::Path;
// use std::path::PathBuf;

use clap::{arg, Arg, Command};
use std::path::PathBuf;

use bio::alphabets::dna;

use std::str::FromStr;

use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;
use noodles::gff::{
    self as gff,
    record::{attributes::Entry, Attributes},
};
use noodles::vcf::{
    self as vcf, header,
    header::record::value::{map::Contig, Map},
    record::reference_bases::Base,
    record::Position,
};
use noodles::vcf::{
    header::info::Key,
    record::chromosome::Chromosome,
    record::info::{field, Field},
    record::Info,
};

use noodles::core::position;
use protein_translate::translate;
use simplelog::*;
use std::cmp::Ordering;

use rayon::prelude::*;
use std::cmp;
use std::fmt;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
pub struct BedEntry {
    // qseqid Query Seq - id
    pub chrom: String,
    // stitle Subject Title
    pub chromStart: u64,
    // sseqid Subject Seq - id
    pub chromEnd: u64,
    // qstart Start of alignment in query*
    pub name: String,
    pub strand: String,
    pub featuretype: tbl::FeatureType,
}

pub struct AllBedEntries {
    pub chrom: String,
    pub entries: Vec<BedEntry>,
}

impl AllBedEntries {
    pub fn join_entries(&mut self, entries_names: Vec<String>) {
        let mut entries = Vec::with_capacity(entries_names.len());

        for locus_tag_name in entries_names.into_iter() {
            let entry_index = self
                .entries
                .iter()
                .position(|x| x.name == locus_tag_name)
                .expect("Could not find the element in cds entries");
            let entry: &BedEntry = self
                .entries
                .get(entry_index)
                .expect("Could not get entry in cds entries");

            entries.push(entry.clone());

            self.entries.remove(entry_index);
        }

        entries.sort_by_key(|r| r.chromStart);

        let first_entry = &entries[0];

        let last_entry = &entries[entries.len() - 1];

        assert_eq!(&first_entry.chrom, &last_entry.chrom);
        // assert_eq!(&first_entry.strand, &last_entry.strand);

        let new_entry = BedEntry {
            chrom: first_entry.chrom.clone(),
            chromStart: first_entry.chromStart,
            chromEnd: last_entry.chromEnd,
            name: first_entry.name.clone(),
            strand: first_entry.strand.clone(),
            featuretype: first_entry.featuretype.clone(),
        };

        self.entries.push(new_entry);
    }

    pub fn sorted(&mut self) {
        self.entries.sort_by_key(|e| e.chromStart)
    }

    pub fn get_by_locus_tag(&self, locus_tag_id: &str) -> Option<&BedEntry> {
        self.entries.iter().find(|&e| e.name == locus_tag_id)
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let seq_len = seq.len();
    let mut vec = Vec::with_capacity(seq_len);

    let mut seqrev = Vec::with_capacity(seq_len);
    for c in seq.iter() {
        seqrev.push(*c);
    }

    seqrev.reverse();

    for nucl in seqrev.iter() {
        let nindex = ASCII_TO_INDEX[*nucl as usize];
        let cbase = NUCL_COMPLEMENT[nindex];
        vec.push(cbase)
    }

    vec
}

#[derive(Debug)]
pub enum StopCodonError {
    PrematureStopCodon(usize),
    WithoutStopCodon,
}

impl std::error::Error for StopCodonError {}

impl fmt::Display for StopCodonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            StopCodonError::PrematureStopCodon(nstop) => {
                write!(
                    f,
                    "CDS with {} premature(s) stop codon(s) within CDS region",
                    nstop
                )
            }
            StopCodonError::WithoutStopCodon => write!(f, "CDS withou a stop codon at its end"),
        }
    }
}

pub fn check_stop_codons(peptide: &str) -> Result<i32, StopCodonError> {
    let n_stop_codons = peptide.as_bytes().iter().filter(|&&x| x == b'*').count();

    if n_stop_codons == 0 {
        return Err(StopCodonError::WithoutStopCodon);
        // eprintln!("WARNING: The created entry cds doesnt have a STOP CODON.");
    } else if n_stop_codons == 1 {
        if *peptide.as_bytes().last().unwrap() == b'*' {
            return Ok(1);
        } else {
            return Err(StopCodonError::WithoutStopCodon);
        }
        // eprintln!("OK: The created entry has only one stop codon.");
    }

    // eprintln!(
    //     "WARNING: The created entry has more than one stop codon (N = {})",
    //     &n_stop_codons
    // );
    Err(StopCodonError::PrematureStopCodon(n_stop_codons))
}

pub fn get_enzyme_map(base_path: &str) -> HashMap<String, String> {
    let mut enzyme_path = PathBuf::from(&base_path);
    // enzyme_path.pop();
    // enzyme_path.push("enzyme.dat");
    let path_enzyme_file = Path::new(&enzyme_path);
    let file_contents = fs::read_to_string(path_enzyme_file).expect("Unable to read file");

    let mut enz_map: HashMap<String, String> = file_contents
        .split("//")
        .map(|snip| {
            let mut ecid = "".to_string();
            let mut ecname = "".to_string();
            for f in snip.lines() {
                if f.starts_with("ID") {
                    ecid = f[5..].to_string();
                }

                if f.starts_with("DE") {
                    ecname.push_str(&f[5..]);
                }
            }
            ecname.pop();
            // let necname = ecname.replace("]", "");
            (ecid, ecname)
        })
        .collect();

    let enz_dict = enz_map.clone();

    for (key, val) in enz_map.iter_mut() {
        let re = Regex::new(r"Transferred entry: (\d+\.\d+\.\d+.\d+)").unwrap();

        if re.find(val).is_some() {
            let caps = re.captures(val);
            let ec_entry = caps.unwrap().get(1).unwrap().as_str();

            *val = enz_dict.get(ec_entry).unwrap().to_string();
        }
    }

    enz_map
}

static NUCL_COMPLEMENT: [u8; 4] = ['T' as u8, 'G' as u8, 'C' as u8, 'A' as u8];

/// Maps an ASCII character to array index
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
static ASCII_TO_INDEX: [usize; 128] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79    (65 = A, 67 = C, 71 = G)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95    (84 = T, 85 = U)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111   (97 = a, 99 = c, 103 = g)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127  (116 = t, 117 = u)
];

pub fn get_genomic_sequence(
    rec_id: &str,
    start: u64,
    end: u64,
    faidx: &mut IndexedReader<File>,
) -> Vec<u8> {
    let mut seq = Vec::new();

    faidx
        .fetch(rec_id, start, end)
        .expect("Couldn't fetch interval");
    faidx.read(&mut seq).expect("Couldn't read the interval");
    seq
}

pub fn get_corrected_tbl(
    res: &mut tbl::Reader,
    hash_tbl_feat: &mut HashMap<(String, tbl::FeatureType), (tbl::Feature, String)>,
    genes: AnnotMap<String, String>,
) -> tbl::Reader {
    for tbl_entry in res.iter_mut() {
        let chr_name = tbl_entry.seqid.clone();
        for tbl_rec in tbl_entry.iter_mut() {
            if tbl_rec.qualifiers.contains_key("locus_tag") {
                hash_tbl_feat.insert(
                    (
                        tbl_rec.qualifiers.get("locus_tag").unwrap().clone(),
                        tbl_rec.feature_key.clone(),
                    ),
                    (tbl_rec.clone(), chr_name.clone()),
                );
            } else {
                // try to catch locus_tag from corresponding gene annotation

                let gene_strand = match tbl_rec.strand {
                    gff::record::Strand::None => bio_types::strand::ReqStrand::Forward,
                    gff::record::Strand::Forward => bio_types::strand::ReqStrand::Forward,
                    gff::record::Strand::Reverse => bio_types::strand::ReqStrand::Reverse,
                    gff::record::Strand::Unknown => bio_types::strand::ReqStrand::Forward,
                };

                let query = bio_types::annot::contig::Contig::new(
                    chr_name.clone(),
                    tbl_rec.start as isize,
                    (tbl_rec.end - tbl_rec.start) as usize,
                    gene_strand,
                );
                let hits: Vec<&String> = genes.find(&query).map(|e| e.data()).collect();
                assert_eq!(hits.len(), 1);
                let locus_tag_found = hits[0].clone();
                tbl_rec
                    .qualifiers
                    .insert("locus_tag".to_string(), locus_tag_found);

                hash_tbl_feat.insert(
                    (
                        tbl_rec.qualifiers.get("locus_tag").unwrap().clone(),
                        tbl_rec.feature_key.clone(),
                    ),
                    (tbl_rec.clone(), chr_name.clone()),
                );

                // panic!("The annotation {:?} doesnt have locus_tag. Please specify locus_tag for all features!", tbl_rec);
            }
        }
    }
    res.clone()
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("liftfix")
        .version("v0.1.0")
        .author("Rodrigo Theodoro Rocha <theodorobiotec@gmail.com>")
        .about(
            "Fix annotations (tbl) based on VCF generated after running frame_stats.
        It uses as input the FASTA produced by bcftools consensus.",
        )
        .arg_required_else_help(true)
        .arg(
            Arg::new("input")
                .value_name("INPUT")
                .help("Input tbl")
                .required(true), // .min_values(1),
        )
        .arg(
            Arg::new("genome")
                .value_name("GENOME")
                .short('r')
                .help("Assembly in FASTA. Requires faidx index.")
                .takes_value(true),
        )
        .arg(
            Arg::new("variants")
                .value_name("VARIANTS")
                .help("VCF (.GZ - FIX THIS).")
                .takes_value(true)
                .short('v')
                .long("variants"),
        )
        .arg(
            Arg::new("ecfile")
                .value_name("ECFILE")
                .help("ec file")
                .takes_value(true)
                .short('e')
                .long("ecfile"),
        )
        .arg(
            Arg::new("chain")
                .value_name("CHAIN")
                .help("Chain file (txt).")
                .takes_value(true)
                .short('c')
                .long("chain"),
        )
        .get_matches();

    let input_path_tbl_str = matches.value_of("input").expect("Specify input argument.");
    let input_path_genome_str = matches
        .value_of("genome")
        .expect("Specify genome argument.");
    let input_path_vcf_str = matches
        .value_of("variants")
        .expect("Specify VCF file in argument.");
    let input_enzyme_str = matches
        .value_of("ecfile")
        .expect("Specify EC file according to EC FTP.");
    let input_path_chain_str = matches
        .value_of("chain")
        .expect("Specify chain file in argument.");

    let input_path = PathBuf::from(&input_path_vcf_str);

    // READ VCF
    let mut vcf_input = File::open(input_path_vcf_str)
        .map(bgzf::Reader::new)
        .map(vcf::Reader::new)?;

    let header: noodles::vcf::header::Header = vcf_input.read_header()?.parse()?;

    // READ GENOME
    let input_path_genome = PathBuf::from(&input_path_genome_str);
    let path_fsa = Path::new(&input_path_genome);
    let mut faidx = IndexedReader::from_file(&path_fsa)?;

    // READ TBL ANNOTATION FILE
    // let input_path_tbl = PathBuf::from(&input_path_tbl_str);
    let mut res = File::open(input_path_tbl_str).map(|mut x| tbl::Reader::new(&mut x))?;
    // let res_first = File::open(input_path_tbl_str).map(|mut x| tbl::Reader::new(&mut x))?;

    let mut hash_tbl_feat = HashMap::new();

    let mut genes: AnnotMap<String, String> = AnnotMap::new();
    for chr_entry in res.iter() {
        for gene_entry in chr_entry
            .iter()
            .filter(|elt| matches!(elt.feature_key, FeatureType::Gene))
        {
            let gene_strand = match gene_entry.strand {
                gff::record::Strand::None => bio_types::strand::ReqStrand::Forward,
                gff::record::Strand::Forward => bio_types::strand::ReqStrand::Forward,
                gff::record::Strand::Reverse => bio_types::strand::ReqStrand::Reverse,
                gff::record::Strand::Unknown => bio_types::strand::ReqStrand::Forward,
            };
            let gene = bio_types::annot::contig::Contig::new(
                chr_entry.seqid.clone(),
                gene_entry.start as isize,
                (gene_entry.end - gene_entry.start) as usize,
                gene_strand,
            );
            genes.insert_at(
                gene_entry.qualifiers.get("locus_tag").unwrap().clone(),
                &gene,
            )
        }
    }

    let res_new = get_corrected_tbl(&mut res, &mut hash_tbl_feat, genes);

    // GET ENTRIES TO JOIN BASED ON VCF ANNOTATION FIELD LTID
    let mut tojoindict = HashMap::new();
    for result in vcf_input.records(&header) {
        let record = result?;
        println!(
            "{:?}",
            record
                .info()
                .get(&Key::Other("LTID".parse().unwrap()))
                .expect("VCF entry doesnt have field which specifies a list of cds to join (LTID)")
                .value()
        );

        let ltid = record
            .info()
            .get(&Key::Other("LTID".parse().unwrap()))
            .expect("VCF entry doesnt have field which specifies a list of cds to join (LTID)")
            .value()
            .expect("Unpack VCF error");

        let lt_join = if let field::Value::String(cds_to_join) = ltid {
            Some(cds_to_join.clone())
        } else {
            None
        };

        if let Chromosome::Name(seqid) = record.chromosome() {
            dbg!(&seqid);
            tojoindict.entry(seqid.clone()).or_insert(Vec::new());
            if let Some(x) = tojoindict.get_mut(&seqid.clone()) {
                x.push(lt_join);
            }
        }
    }

    let mut chain_file = File::open(input_path_chain_str)?;
    let mut chainout = String::new();
    chain_file.read_to_string(&mut chainout)?;

    dbg!("START CHAIN PARSE FILE");

    let (_, res_chain) = chain::parse_file(&chainout).expect("Couldn't parse chain file.");

    dbg!("END CHAIN PARSE FILE");

    let chains: AnnotMap<String, LiftBlockBuild> = res_chain.to_lift();

    dbg!("CHAINED TO LIFT");

    let mut bedouts: Vec<AllBedEntries> = Vec::new();
    for tbl_entry in res_new.iter() {
        dbg!("BED IN = ", &tbl_entry);
        let mut outbed = Vec::new();
        for tbl_rec in tbl_entry.iter() {
            let contig = tbl_entry.seqid.clone();
            let start_coordinate = tbl_rec.start - 1;
            let end_coordinate = tbl_rec.end;

            let query = Pos::new(
                tbl_entry.seqid.clone(),
                start_coordinate as isize,
                ReqStrand::Forward,
            );

            let hits = chains
                .find(&query)
                .map(|e| e.data())
                .collect::<Vec<&LiftBlockBuild>>();

            // println!("BLOCK LEN = {:?}", &hits[0].alignment_op.len());

            dbg!("LIFTING?");
            // dbg!(chains.len());
            let lifted_pos_start = hits[0].liftover2assembly(start_coordinate);
            let lifted_pos_end = hits[0].liftover2assembly(end_coordinate);

            let feat_strand = match &tbl_rec.strand {
                noodles::gff::record::Strand::None => ".".to_string(),
                noodles::gff::record::Strand::Forward => "+".to_string(),
                noodles::gff::record::Strand::Reverse => "-".to_string(),
                noodles::gff::record::Strand::Unknown => ".".to_string(),
            };

            // println!("{:?}", &tbl_rec);
            // println!(
            //     "{}\t{}\t{}\t{}\t.\t{}",
            //     &contig,
            //     &lifted_pos_start.unwrap(),
            //     &lifted_pos_end.unwrap(),
            //     &tbl_rec.qualifiers.get("locus_tag").unwrap(),
            //     &feat_strand
            // );

            let cds_bed_entry = BedEntry {
                chrom: contig,
                chromStart: lifted_pos_start.unwrap(),
                chromEnd: lifted_pos_end.unwrap(),
                name: tbl_rec.qualifiers.get("locus_tag").unwrap().clone(),
                strand: feat_strand,
                featuretype: tbl_rec.feature_key.clone(),
            };

            outbed.push(cds_bed_entry);
        }

        let mut bedout = AllBedEntries {
            chrom: outbed[0].chrom.clone(),
            entries: outbed,
        };
        let tojoinentries = tojoindict
            .get(&bedout.chrom)
            .expect("Could not find entry among entries to join");

        let mut entries_names_to_join = vec![];
        let mut toavoiddup = vec![];
        for tj in tojoinentries.iter() {
            let tjc = tj.clone().expect("ltid entry");
            if !toavoiddup.iter().any(|e| e == &tjc) {
                let vec_cds_to_join: Vec<String> = tjc.split(',').map(|x| x.to_string()).collect();
                entries_names_to_join.push(vec_cds_to_join);
                toavoiddup.push(tjc.clone());
            }
        }

        for e in entries_names_to_join.into_iter() {
            dbg!(&e);
            let joined_names = &e.join(",");
            let out_locus_tag = &e[0].clone();
            bedout.join_entries(e);

            let new_entry = bedout
                .get_by_locus_tag(out_locus_tag)
                .expect("Couldnt find newly created joined entry");
            eprintln!(
                "Joining CDS entries {} into entry within locus_tag={}",
                &joined_names, &out_locus_tag
            );

            let cds_genomic_seq = get_genomic_sequence(
                &new_entry.chrom,
                new_entry.chromStart,
                new_entry.chromEnd,
                &mut faidx,
            );

            let peptide = if new_entry.strand == "-" {
                let rc = reverse_complement(&cds_genomic_seq);
                translate(&rc)
            } else {
                translate(&cds_genomic_seq)
            };

            let n_stop_codons = peptide.as_bytes().iter().filter(|&&x| x == b'*').count();

            dbg!(peptide);
            dbg!(new_entry);
            dbg!(&n_stop_codons);

            if n_stop_codons == 0 {
                eprintln!("WARNING: The created entry cds doesnt have a STOP CODON.");
            } else if n_stop_codons == 1 {
                eprintln!("OK: The created entry has only one stop codon.");
            } else {
                eprintln!(
                    "WARNING: The created entry has more than one stop codon (N = {})",
                    &n_stop_codons
                );
            }
        }

        bedout.sorted();
        bedouts.push(bedout);
    }

    for c in bedouts.iter() {
        for cc in c.entries.iter() {
            let cds_genomic_seq =
                get_genomic_sequence(&cc.chrom, cc.chromStart, cc.chromEnd, &mut faidx);

            let peptide = if cc.strand == "-" {
                let rc = reverse_complement(&cds_genomic_seq);
                translate(&rc)
            } else {
                translate(&cds_genomic_seq)
            };

            match check_stop_codons(&peptide) {
                Ok(_) => {
                    dbg!("TUDO OK", &peptide);
                    ()
                }
                Err(e) => match e {
                    StopCodonError::PrematureStopCodon(_) => {
                        dbg!("PREMATURE", &peptide);
                    }
                    StopCodonError::WithoutStopCodon => {
                        dbg!("SEM STOP", &peptide);
                    }
                },
            }
            // dbg!(peptide);
            // dbg!(cc);
        }
    }

    let mut output_path = PathBuf::from(&input_path_tbl_str);
    output_path.set_extension("gff3");

    let mut writer = gff::Writer::new(Vec::new());
    let version = gff::Directive::GffVersion(Default::default());
    writer.write_directive(&version)?;

    let enz_map = get_enzyme_map(&input_enzyme_str);

    for chrom_instance in bedouts.iter() {
        for bed_entry in chrom_instance.entries.iter() {
            let feat = &hash_tbl_feat
                .get(&(bed_entry.name.clone(), bed_entry.featuretype.clone()))
                .expect("Could not find TBL entry in hashmap")
                .0;
            let mut vec_attributes: Vec<Entry> = Vec::new();
            // dbg!(&feat.qualifiers);
            let mut new_attributes_dict = feat.qualifiers.clone();

            for (k, v) in &feat.qualifiers {
                vec_attributes.push(Entry::new(k.clone(), v.clone()));
            }

            // Correct EC NUMBERS according to ENZYME DB
            if feat.qualifiers.contains_key("EC_number") {
                let ec_value = feat.qualifiers.get("EC_number");
                let ec_enzymap_name = match ec_value {
                    Some(ec_number) => enz_map.get(ec_number),
                    None => None,
                };

                // If enzyme ec number couldnt be found in ENZYME MAP DICT, thus
                // lets remove it from final annotation attributes
                if ec_enzymap_name.is_none() {
                    new_attributes_dict.remove("EC_number");
                    let mut ec_value_raw: String = String::new();
                    if let Some(ecv) = ec_value {
                        ec_value_raw.push_str(ecv);
                    }
                    let mut value_to_note = "EC_number_".to_string();
                    value_to_note.push_str(&ec_value_raw);
                    new_attributes_dict.append("Note".to_string(), value_to_note.to_string());
                } else {
                    let val = new_attributes_dict.get_mut("EC_number").unwrap().clone();

                    let _ = new_attributes_dict.remove_entry("EC_number");
                    let _ = new_attributes_dict.remove_entry("product");
                    if let Some(product_name) = enz_map.get(&val) {
                        new_attributes_dict.insert("product".to_string(), product_name.to_string());
                    }
                    new_attributes_dict.insert("ec_number".to_string(), val);
                }
            }

            if matches!(bed_entry.featuretype, FeatureType::CDS) {
                // Correct PSEUDOGENES according to NCBI submission guidelines
                // https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
                let cds_genomic_seq = get_genomic_sequence(
                    &bed_entry.chrom,
                    bed_entry.chromStart,
                    bed_entry.chromEnd,
                    &mut faidx,
                );

                let peptide = if bed_entry.strand == "-" {
                    let rc = reverse_complement(&cds_genomic_seq);
                    translate(&rc)
                } else {
                    translate(&cds_genomic_seq)
                };

                match check_stop_codons(&peptide) {
                    Ok(_) => (),
                    Err(e) => match e {
                        StopCodonError::PrematureStopCodon(_) => {
                            dbg!("PREMATURE", &peptide);

                            new_attributes_dict.insert("pseudo".to_string(), "true".to_string());
                        }
                        StopCodonError::WithoutStopCodon => {
                            dbg!("SEM STOP", &peptide);

                            new_attributes_dict.insert("pseudo".to_string(), "true".to_string());
                        }
                    },
                }
            }

            // Transform all the attributes that are not listed as
            // - locus_tag
            // - transcript_id
            // - protein_id
            // - product
            // - pseudo
            // - pseudogene
            // - DBxref
            // - ec_number
            // - gene
            // - gene_synonym
            // - description
            // - exception
            // - transl_except
            // - function
            // - experiment
            // - old_locus_tag
            // - mobile_element
            // - ncRNA_class
            // - regulatory_class
            // - recombination_class
            // To Notes attribute

            // TO DO!!!

            // Initiate GFF entry records build from updated atrributes
            // that were corrected
            // they are: EC number failures, stop codon inside genes
            let mut vec_attributes: Vec<Entry> = Vec::new();
            // dbg!(&feat.qualifiers);

            for (k, v) in &new_attributes_dict {
                vec_attributes.push(Entry::new(k.clone(), v.clone()));
            }

            let mut gff_record: gff::Record;

            gff_record = gff::Record::builder()
                .set_reference_sequence_name(bed_entry.chrom.clone())
                .set_start(position::Position::new((bed_entry.chromStart + 1) as usize).unwrap())
                .set_end(position::Position::new(bed_entry.chromEnd as usize).unwrap())
                .set_type(feat.feature_key.to_string())
                .set_strand(noodles::gff::record::Strand::from_str(&bed_entry.strand).unwrap())
                .set_attributes(Attributes::from(vec_attributes))
                .build();

            writer
                .write_record(&gff_record)
                .expect("Error writing gff record");
            println! {"{:?}\n\n", &feat}
        }
    }

    let mut file = File::create(output_path)?;
    // let writer_file = BufWriter::new(file);
    file.write_all(writer.get_ref())
        .expect("Unable to write data.");

    // println!("{:?}", &chains);

    // dbg!(enz_map);

    Ok(())
}
