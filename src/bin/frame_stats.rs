/*!

# What is it

Search for possible frameshift errors in coding regions.

### Usage, option summary and outputs

```text
USAGE:
    frame_stats [OPTIONS] <INPUT>

ARGS:
    <INPUT>    Input tbl

OPTIONS:
    -d, --diamond <DIAMOND>    Diamond result after running blastx with gene sequences
    -g, --genes <GENES>        Fasta file with genes sequences (DNA). The gene id must match,
                               locus_tag in tbl input.
    -h, --help                 Print help information
    -r <GENOME>                Assembly in FASTA. Requires faidx index.
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

use bio::io::fasta::{self, IndexedReader, Record};
use noodles::bgzf::{self as bgzf};
use std::io;
use std::io::prelude::*;
// use nom::error::dbg_dmp;
// use rs_genomesub::tbl::Feature;
// use uniprot::uniprot::DbReference;
use ureq::Response;

extern crate rs_genomesub;
use crate::rs_genomesub::blast_tabular::*;
use crate::rs_genomesub::tbl;

use itertools::Itertools;
use std::collections::{BTreeMap, HashMap};

use std::fs;
use std::fs::File;
// use std::io::Write;
use std::path::Path;
// use std::path::PathBuf;

use clap::{arg, Arg, Command};
use std::path::PathBuf;

use bio::alphabets::dna;

use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;
use noodles::vcf::{
    self as vcf, header,
    header::record::value::{map::Contig, Map},
    record::reference_bases::Base,
    record::Position,
};
use noodles::vcf::{
    header::info::Key,
    record::info::{field, Field},
    record::Info,
};
use protein_translate::translate;
use simplelog::*;
use std::cmp::Ordering;

use rayon::prelude::*;
use std::cmp;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};

/// Given a Protein ID, get the corresponding genomic sequence through an API request to ENA server.
pub fn get_record_from_ena(response: Response) -> Result<Record, std::io::Error> {
    let reader = libflate::gzip::Decoder::new(response.into_reader()).unwrap();

    let entry = uniprot::uniprot::parse(std::io::BufReader::new(reader)).next();

    // CHECK ABOVE

    if entry.is_some() {
        let mut prop_ena_id: Option<String> = None;
        let entry_res = entry.unwrap();
        if entry_res.is_ok() {
            for db_ref in entry_res.unwrap().db_references.iter() {
                for prop in db_ref.property.iter() {
                    if prop.value == "Genomic_DNA" {
                        // dbg!(&db_ref);
                        for prop_id in db_ref.property.clone().iter() {
                            if prop_id.ty == "protein sequence ID" {
                                prop_ena_id = Some(prop_id.value.clone());
                            }
                        }
                    }
                }
            }

            let query_url_ena = format!(
                "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/ena_sequence/{}/fasta",
                prop_ena_id.unwrap()
            );
            let req_ena = ureq::get(&query_url_ena);

            let reader_ena = fasta::Reader::new(req_ena.call().unwrap().into_reader());

            return reader_ena.records().next().unwrap();
        }
    }

    Err(std::io::Error::new(std::io::ErrorKind::Other, "server"))
}

pub fn semiglobal_align_both_directions(
    boundary_seq: &Vec<u8>,
    ena_seq: &Vec<u8>,
) -> Vec<(i32, Alignment, f64)> {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

    let alignments = (0..=1)
        .into_par_iter()
        .map(|strand| {
            let target_seq = if strand == 0 {
                ena_seq.clone()
            } else {
                dna::revcomp(ena_seq)
            };

            let mut aligner =
                Aligner::with_capacity(boundary_seq.len(), ena_seq.len(), -5, -1, &score);

            let alignment = aligner.semiglobal(boundary_seq, &target_seq);

            let (match_count, mismatch_count) =
                alignment
                    .operations
                    .iter()
                    .fold((0, 0), |(matches, mismatches), x| match x {
                        Match => (matches + 1, mismatches),
                        Subst => (matches, mismatches + 1),
                        Del => (matches, mismatches),
                        Ins => (matches, mismatches),
                        Xclip(_) => (matches, mismatches),
                        Yclip(_) => (matches, mismatches),
                    });

            let gap_free_score =
                (match_count as f64) / (match_count as f64 + mismatch_count as f64);

            (strand, alignment, gap_free_score)
        })
        .collect::<Vec<(i32, Alignment, f64)>>();

    alignments
}

// alt in ["A", "C", "G", "T", "-"]
pub fn expand_homopolymer(pos: usize, alt: &[u8], reference: &[u8]) -> usize {
    let mut i_threeprime = 0;
    // let mut repeated_bases = Vec::new();
    // i put a pos
    while pos + i_threeprime < reference.len() && ([reference[pos + i_threeprime]] == *alt) {
        i_threeprime += 1;
    }

    let mut i_fiveprime = 0;
    // let mut repeated_bases = Vec::new();
    let mut ref_pos = pos.checked_sub(i_fiveprime + 1);
    while let Some(rp) = ref_pos {
        if ([reference[rp]] == *alt) {
            i_fiveprime += 1;
            ref_pos = pos.checked_sub(i_fiveprime + 1);
        } else {
            ref_pos = None;
        }
    }

    i_threeprime + i_fiveprime
}

// pub fn get_record_from_insilico_translation(
//     response: Response,
//     prot_id: &str,
// ) -> Result<Record, std::io::Error> {
//     let reader = libflate::gzip::Decoder::new(response.into_reader()).unwrap();

//     let entry = uniprot::uniref::parse(std::io::BufReader::new(reader))
//         .next()
//         .unwrap();

//     let protseq = if let uniprot::uniref::Member { sequence: seq, .. } =
//         entry.unwrap().representative_member
//     {
//         if let Some(uniprot::uniparc::Sequence {
//             sequence: string_sequence,
//             ..
//         }) = seq
//         {
//             string_sequence
//         } else {
//             unimplemented!()
//         }
//     } else {
//         unimplemented!()
//     };

//     println!("STRING SEQ: {}", &protseq);

//     // println!("{:?}", &entry);

//     // let protseq = b"ATC";
//     // let read_id = "ehissoai";
//     let read_id = format!("InSilico_{}", prot_id);
//     let protseq_u8 = protseq.as_bytes();
//     let dna = translate(protseq_u8);

//     dbg!("UAI, ", &dna);

//     let description = None;
//     let record = Record::with_attrs(&read_id, description, dna.as_bytes());

//     dbg!("{:?}", &record);

//     Ok(record)
// }

pub fn get_genomic_sequence(
    rec_id: &str,
    start: u64,
    end: u64,
    faidx: &mut IndexedReader<File>,
) -> Vec<u8> {
    let mut seq = Vec::new();

    faidx
        .fetch(rec_id, start - 1, end)
        .expect("Couldn't fetch interval");
    faidx.read(&mut seq).expect("Coldn't read the interval");
    seq
}

fn main() {
    // let cds_fasta = include_str!("../../genes_out.fa");

    // let input_path = PathBuf::from(&cds_fasta);

    let matches = Command::new("tbl_frames")
        .version("v0.1.0")
        .author("Rodrigo Theodoro Rocha <theodorobiotec@gmail.com>")
        .about("Inspect candidate frame shifts erros within coding regions")
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
            Arg::new("genes")
                .value_name("GENES")
                .help("Fast file with genes sequences (DNA). The gene id must match, locus_tag in tbl input.")
                .takes_value(true)
                .short('g')
                .long("genes"),
        )
        .arg(
            Arg::new("diamond")
                .value_name("DIAMOND")
                .help("Diamond result after running blastx with gene sequences")
                .takes_value(true)
                .short('d')
                .long("diamond"),
        )
        .get_matches();

    let config = ConfigBuilder::new()
        .set_level_color(Level::Error, Some(Color::Rgb(191, 0, 0)))
        .set_level_color(Level::Warn, Some(Color::Rgb(255, 127, 0)))
        .set_level_color(Level::Info, Some(Color::Rgb(192, 192, 0)))
        .set_level_color(Level::Debug, Some(Color::Rgb(63, 127, 0)))
        .set_level_color(Level::Trace, Some(Color::Rgb(127, 127, 255)))
        .set_time_to_local(true)
        .set_time_format("%Y-%m-%d %H:%M:%S".to_string())
        .set_thread_level(LevelFilter::Info)
        .set_thread_mode(ThreadLogMode::Both)
        .build();

    let input_path_tbl_str = matches.value_of("input").unwrap();
    let input_path_genome_str = matches.value_of("genome").unwrap();
    let input_path_genes_str = matches.value_of("genes").unwrap();
    let input_path_diamond_str = matches.value_of("diamond").unwrap();

    let input_path = PathBuf::from(&input_path_genes_str);

    let mut input = File::open(input_path).unwrap();
    let res = fasta::Reader::new(&mut input);

    let out_vcf = File::create("out.vcf.gz")
        .map(bgzf::Writer::new)
        .expect("Unable to create file to write output sequences.");

    let mut out_log = PathBuf::from(&input_path_genes_str);
    out_log.set_file_name("framerust.log");
    let out_log_file = Path::new(&out_log);

    CombinedLogger::init(vec![
        TermLogger::new(
            LevelFilter::Info,
            config,
            TerminalMode::Stdout,
            ColorChoice::Auto,
        ),
        WriteLogger::new(
            LevelFilter::Info,
            Config::default(),
            File::create(out_log_file).unwrap(),
        ),
    ])
    .unwrap();

    let mut writer_vcf = vcf::Writer::new(out_vcf);

    let mut hash_cds_fasta = HashMap::new();

    let mut record;
    let mut record_id;
    for s in res.records() {
        record = s.unwrap();
        record_id = record.id().to_string();
        hash_cds_fasta.insert(record_id, record);
    }

    let input_path_tbl = PathBuf::from(&input_path_tbl_str);
    let mut input_tbl = File::open(input_path_tbl).unwrap();
    let res = tbl::Reader::new(&mut input_tbl);
    let mut hash_tbl_feat = HashMap::new();

    for tbl_entry in res.iter() {
        for tbl_rec in tbl_entry.iter() {
            if tbl_rec.qualifiers.contains_key("locus_tag") {
                hash_tbl_feat.insert(
                    tbl_rec.qualifiers.get("locus_tag").unwrap(),
                    (tbl_rec, tbl_entry.seqid.clone()),
                );
            }
        }
    }

    info!(
        "{} genes sequences (CDS) read from {}.",
        hash_cds_fasta.len(),
        &input_path_genes_str
    );

    // let mut output_fasta: Vec<fasta::Record> = Vec::new();
    let mut output_fasta_ref = Arc::new(Mutex::new(Vec::new()));
    let mut vcf_records_to_write_ref = Arc::new(Mutex::new(Vec::new()));
    let mut entries_no_dna_ref = Arc::new(Mutex::new(Vec::new()));
    let mut complex_loci_ref = Arc::new(Mutex::new(0));
    let mut hash_alignment_pretty_ref =
        Arc::new(Mutex::new(BTreeMap::<(String, usize), String>::new()));
    let mut pair_with_indel_ref = Arc::new(Mutex::new(0));

    let input_path_diamond = PathBuf::from(&input_path_diamond_str);
    let fblastout = fs::read_to_string(input_path_diamond).expect("Unable to read diamond result");
    // let fblastout = include_str!("../../examples/output_diamond_sample.txt");

    // Parse diamond result into BlastFeature struct
    let (_inp, res) = parse_input(&fblastout).unwrap();

    let input_path_genome = PathBuf::from(&input_path_genome_str);
    let path_fsa = Path::new(&input_path_genome);

    let mut faidx = IndexedReader::from_file(&path_fsa).unwrap();

    std::env::set_var("RAYON_NUM_THREADS", "8");

    // let count_possible_frameshifts = res
    //     .into_iter()
    //     .tuple_windows()
    //     .filter(|(s1, s2)| s1.sseqid == s2.sseqid)
    //     .count();
    let resv2 = res.clone();
    let count_possible_frameshifts = res
        .into_iter()
        .tuple_windows()
        .filter(|(s1, s2)| (s1.sseqid == s2.sseqid))
        .collect::<Vec<(BlastFeature, BlastFeature)>>();

    info!("Starting <green>{} workers</>.", 8);
    info!("Detecting adjacent genes on the same strand that give hits against the same subject (common BLAST hit).");
    let mut count_possible_frameshifts_v2 = resv2
        .iter()
        .group_by(|x| (x.sseqid.clone()))
        .into_iter()
        .map(|(sseqid, group)| MultipleBlastFeatures {
            id: sseqid,
            bfeatures: group.map(|u| u.clone()).collect(),
        })
        .filter(|s| s.bfeatures.len() > 1)
        .collect::<Vec<MultipleBlastFeatures>>();

    let count_possible_frameshifts_v3 = resv2
        .iter()
        .group_by(|x| (x.sseqid.clone()))
        .into_iter()
        .map(|(sseqid, group)| MultipleBlastFeatures {
            id: sseqid,
            bfeatures: group.map(|u| u.clone()).collect(),
        })
        .filter(|s| {
            if s.bfeatures.len() == 1 {
                let qlen = (s.bfeatures[0].qlen / 3) as f64;
                let slen = s.bfeatures[0].slen as f64;

                if qlen / slen < 0.95 {
                    true
                } else {
                    false
                }
            } else {
                false
            }
        })
        .collect::<Vec<MultipleBlastFeatures>>();

    let count_possible_frameshifts_num = &count_possible_frameshifts.len();
    eprintln!(
        "Quantity of pairs of annotated coding sequences (CDS) matching the same uniref entry: {}.\nThis fact suggests that there may be frameshifts interrupting a pair.",
        &count_possible_frameshifts_num
    );

    let count_possible_frameshifts_v2_num = &count_possible_frameshifts_v2.len();
    eprintln!(
        "V2: Quantity of pairs of annotated (AGGREGATED) coding sequences (CDS) matching the same uniref entry: {}.\nThis fact suggests that there may be frameshifts interrupting a pair.",
        &count_possible_frameshifts_v2_num
    );

    let count_possible_frameshifts_v3_num = &count_possible_frameshifts_v3.len();
    eprintln!(
        "V3: Quantity of single features wich has qlen/slen < 0.95: {}.",
        &count_possible_frameshifts_v3_num
    );

    // let mut boundary_seq_description_v2 = format!("{}:{}-{}", chr_id, tbl_bf1.start, tbl_bf2.end);

    count_possible_frameshifts_v2
        .par_iter_mut()
        .for_each(|pfs| {
            let mut faidx_child = IndexedReader::from_file(&path_fsa).unwrap();
            let output_fasta_child = output_fasta_ref.clone();
            let vcf_records_to_write_child = vcf_records_to_write_ref.clone();
            let entries_no_dna_child = entries_no_dna_ref.clone();
            let complex_loci_child = complex_loci_ref.clone();
            let hash_alignment_pretty_child = hash_alignment_pretty_ref.clone();
            let pair_with_indel_child = pair_with_indel_ref.clone();

            pfs.bfeatures.sort_by(|a, b| {
                let a_start = hash_tbl_feat[&a.qseqid].0.start;
                let b_start = hash_tbl_feat[&b.qseqid].0.start;
                if (a_start < b_start) {
                    Ordering::Less
                } else if (a_start == b_start) {
                    Ordering::Equal
                } else {
                    Ordering::Greater
                }
            });

            let bf_names: Vec<String> = pfs.bfeatures.iter().map(|f| f.qseqid.clone()).collect();

            let first_bf_feature = pfs.bfeatures.first().unwrap();
            let last_bf_feature = pfs.bfeatures.last().unwrap();

            let (_, chr_id) = hash_tbl_feat.get(&first_bf_feature.qseqid).unwrap();
            let genome_base_start = hash_tbl_feat[&first_bf_feature.qseqid].0.start;
            let genome_base_end = hash_tbl_feat[&last_bf_feature.qseqid].0.end;

            let boundary_seq =
                get_genomic_sequence(chr_id, genome_base_start, genome_base_end, &mut faidx_child);

            let read_id = format!("my_genomic_boundary_{}", bf_names.join("_"));
            let boundary_seq_description = format!("{}:{}", chr_id, bf_names.join("-"));
            let description = Some(&boundary_seq_description[..]);
            // let sequence = b"ACGT";
            let record = Record::with_attrs(&read_id[..], description, &boundary_seq);

            if let Some(query) = first_bf_feature.sseqid.strip_prefix("UniRef90_") {
                let query_url = format!(
                    "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
                    // "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
                    &query
                );

                println!("{}", &query_url);

                let visual_url = format!("https://www.uniprot.org/uniprot/{}", &query);

                debug!(
                "Performing an HTTP GET request against {} to retrieve DNA sequence for entry {}.",
                &visual_url, &first_bf_feature.sseqid
                );
                let req = ureq::get(&query_url)
                    .set("Accept", "application/xml")
                    .call();

                let reader = match req {
                    Ok(r) => get_record_from_ena(r),
                    Err(e) => {
                        let query_url = format!(
                        // "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
                        "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
                        query
                    );

                        let req = ureq::get(&query_url)
                            .set("Accept", "application/xml")
                            .call()
                            .unwrap();

                        warn!(
                            "Could not retrieve DNA sequence from entry UniRef90_{}",
                            &query
                        );
                        // get_record_from_insilico_translation(req, query)
                        // CHECK
                        // unimplemented!()
                        Err(std::io::Error::new(
                            std::io::ErrorKind::Other,
                            "unretrieved",
                        ))
                    }
                };

                // let ena_rec = match reader {
                //     Ok(rec) => {
                //         // output_fasta.push(rec.clone());
                //         Ok(rec)
                //     }
                //     Err(e) => {
                //         eprintln!("Could not find Record");
                //         Err(e)
                //     }
                // };

                if reader.is_ok() {
                    let ena = reader.unwrap();

                    let ena_seq = ena.seq().to_vec();
                    let alignments = semiglobal_align_both_directions(&boundary_seq, &ena_seq);

                    // Returns the element that gives the maximum gap_free_score value among forward and reverse alignments.
                    let best_aln = alignments
                        .iter()
                        .max_by(|(_, _, a_gap_free_score), (_, _, b_gap_free_score)| {
                            a_gap_free_score.partial_cmp(b_gap_free_score).unwrap()
                        })
                        .unwrap();

                    let pos_indels = best_aln
                        .1
                        .operations
                        .iter()
                        .enumerate()
                        .filter(|(_, &s)| match s {
                            Ins => true,
                            Del => true,
                            _ => false,
                        })
                        .map(|(idx, op)| {
                            let best_aln_cp = best_aln.1.clone();
                            let dels_in_opvec = &best_aln_cp.operations[..idx]
                                .into_iter()
                                .filter(|el| match el {
                                    Del => true,
                                    _ => false,
                                })
                                .count();

                            let ins_in_opvec = &best_aln_cp.operations[..idx]
                                .into_iter()
                                .filter(|el| match el {
                                    Ins => true,
                                    _ => false,
                                })
                                .count();

                            let var_pos = (idx + genome_base_start as usize) - dels_in_opvec;

                            let strand = best_aln.0;
                            // extract indel first base
                            // strand == 0 -> Forward
                            let var_char = match op {
                                Ins => {
                                    let var_pos =
                                        (idx + genome_base_start as usize) - dels_in_opvec;

                                    let var_pos_0 = var_pos;
                                    let var_pos_1 = var_pos + 1;
                                    let variant_char = if strand == 0 {
                                        get_genomic_sequence(
                                            chr_id,
                                            (var_pos_0 + 1) as u64,
                                            (var_pos_1) as u64,
                                            &mut faidx_child,
                                        )
                                    } else {
                                        get_genomic_sequence(
                                            chr_id,
                                            (var_pos_0) as u64,
                                            (var_pos_1 - 1) as u64,
                                            &mut faidx_child,
                                        )
                                    };
                                    variant_char
                                }
                                Del => {
                                    let var_pos =
                                        (idx + genome_base_start as usize) - dels_in_opvec;

                                    let var_pos_0 = var_pos;
                                    let var_pos_1 = var_pos + 1;
                                    let variant_char_old = if strand == 0 {
                                        get_genomic_sequence(
                                            chr_id,
                                            (var_pos_0 + 1) as u64,
                                            (var_pos_1) as u64,
                                            &mut faidx_child,
                                        )
                                    } else {
                                        get_genomic_sequence(
                                            chr_id,
                                            (var_pos_0) as u64,
                                            (var_pos_1 - 1) as u64,
                                            &mut faidx_child,
                                        )
                                    };
                                    let variant_char = get_genomic_sequence(
                                        chr_id,
                                        (var_pos_0) as u64,
                                        (var_pos_1 - 1) as u64,
                                        &mut faidx_child,
                                    );
                                    variant_char
                                }
                                _ => panic!("Must rely on INDELS"),
                            };

                            let var_pos_bhit = match op {
                                Ins => var_char[0],
                                Del => {
                                    let var_pos_bhit = idx - ins_in_opvec;
                                    let variant_char_bhit_in = if strand == 0 {
                                        ena_seq[var_pos_bhit]
                                    } else {
                                        let ena_seq_rc = dna::revcomp(&ena_seq);
                                        ena_seq_rc[var_pos_bhit]
                                    };
                                    variant_char_bhit_in
                                }
                                _ => panic!("Must rely on INDELS"),
                            };

                            (op, var_pos, var_pos + 1, var_char[0], idx, var_pos_bhit)
                        })
                        .collect::<Vec<(&AlignmentOperation, usize, usize, u8, usize, u8)>>();

                    let mut data_grouped = Vec::new();
                    let mut hs_rl = HashSet::new();
                    let all = HashSet::from_iter(pos_indels.iter().cloned());
                    let mut pos_indels_sorted: Vec<(
                        &AlignmentOperation,
                        usize,
                        usize,
                        u8,
                        usize,
                        u8,
                    )> = pos_indels.iter().copied().collect();
                    pos_indels_sorted.sort_by(|a, b| a.1.cmp(&b.1));
                    // dbg!(&pfs.id, &pos_indels_sorted);
                    for (key, group) in &pos_indels_sorted
                        .clone()
                        .into_iter()
                        .tuple_windows::<(_, _)>()
                        .group_by(|(elt_before, elt_after)| elt_after.4 == (1 + elt_before.4))
                    {
                        // dbg!("KEEEEY?", &pfs.id);
                        let mut hs_inner = HashSet::new();
                        let gc = group.collect::<Vec<(_, _)>>();
                        if key {
                            for g in gc.iter() {
                                // dbg!(&g);
                                let (a, b) = g;

                                hs_rl.insert(a.clone());
                                hs_rl.insert(b.clone());

                                hs_inner.insert(a.clone());
                                hs_inner.insert(b.clone());
                            }
                            // evens stores contiguous runs of indels (i.e, a stretch)
                            let mut evens = hs_inner.into_iter().collect::<Vec<_>>();
                            evens.sort_by(|(_, a, _, _, _, _), (_, b, _, _, _, _)| a.cmp(b));

                            if evens.len() <= 3 {
                                data_grouped.push((key, evens));
                            }
                        }
                    }

                    let mut singletons: HashSet<_> = all
                        .difference(&hs_rl)
                        .map(|&x| {
                            let strand = best_aln.0;
                            // extract indel first base
                            // strand == 0 -> Forward

                            (x.0, x.1, x.2, x.3, x.4, x.5)
                        })
                        .collect();
                    // dbg!(&singletons);

                    let strand = best_aln.0;

                    let mut complex_loci = complex_loci_child.lock().expect("Error on fasta");

                    if !data_grouped.is_empty() {
                        *complex_loci += 1;

                        dbg!("IS HERE");
                        // dbg!(&pfs.id, &singletons);
                        dbg!(&pfs.id, &data_grouped);

                        // code to check if all nucleotides found in a run length are equal
                        // this object only keeps the indels that have same nucleotide
                        // it can be thought as a way to ignore larger indels that are not homopolymers in full length
                        let mut db_nucl: Vec<Option<(_, _, _, _, _, _)>> = data_grouped
                            .iter()
                            .map(|(key, evens)| {
                                // we can do unwrap because we've already checked that the set is not empty
                                // therefore it has at least one element
                                let mut evens_iter = evens.into_iter();
                                let first_el = evens_iter.next().unwrap();
                                let first_nucl = first_el.5;

                                let alignment_op = first_el.0;
                                let idx_aln = first_el.4;

                                evens_iter
                                    .all(|elem| {
                                        let elem_nucl = elem.5;

                                        elem_nucl == first_nucl
                                    })
                                    .then(|| {
                                        let start_position = evens
                                            .iter()
                                            .enumerate()
                                            .map(|(i, (_, start_base, _, _, _, _))| start_base + i)
                                            .min()
                                            .unwrap();
                                        let end_position = evens
                                            .iter()
                                            .enumerate()
                                            .map(|(i, (_, _, end_base, _, _, _))| end_base + i)
                                            .max()
                                            .unwrap();
                                        (
                                            alignment_op,
                                            start_position.clone(),
                                            end_position.clone(),
                                            first_nucl,
                                            idx_aln,
                                            first_nucl,
                                        )
                                    })
                            })
                            .collect();

                        // let db_nucl_count: Vec<Vec<(usize, u8)>> = data_grouped
                        //     .iter()
                        //     .map(|(key, evens)| {
                        //         // we can do unwrap because we've already checked that the set is not empty
                        //         // therefore it has at least one element
                        //         let mut evens_iter = evens.into_iter();

                        //         let count: Vec<(usize, u8)> = evens_iter
                        //             .dedup_by_with_count(|x, y| x.3 == y.3)
                        //             .map(|(c, el)| (c, el.3))
                        //             .collect();

                        //         count
                        //     })
                        //     .collect();

                        // dbg!(&db_nucl_count);

                        for runlength_var in db_nucl.iter() {
                            if runlength_var.is_some() {
                                let rl = runlength_var.clone().unwrap();
                                singletons.insert(rl);
                            }
                        }

                        dbg!(&db_nucl);
                    }

                    // dbg!(&pos_indels);

                    // let mut singletons_vec: Vec<_> = singletons.into_iter().collect();

                    let mut pair_with_indel = pair_with_indel_child.lock().expect("Error on fasta");
                    if !singletons.is_empty() {
                        *pair_with_indel += 1;
                    }

                    // Iterate over each indel observed in best alignment
                    for indel in singletons.iter() {
                        // Indel start position
                        let pos_var = (indel.1 + 1) as usize;
                        let var_idx = indel.4;
                        let alignment_line_pos = var_idx.div_euclid(100);

                        let ena_seq_to_print = if strand == 0 {
                            ena_seq.clone()
                        } else {
                            dna::revcomp(&ena_seq)
                        };

                        let line_align: Vec<String> = best_aln
                            .1
                            .pretty(&boundary_seq, &ena_seq_to_print)
                            .split("\n\n\n")
                            .collect::<Vec<&str>>()
                            .iter()
                            .enumerate()
                            .filter_map(|(idx, l)| {
                                if idx == alignment_line_pos {
                                    Some(l.to_string())
                                } else {
                                    None
                                }
                            })
                            .collect();

                        let mut output_alignment_line = hash_alignment_pretty_child
                            .lock()
                            .expect("Error on writing print alignment");
                        output_alignment_line
                            .insert((chr_id.clone(), pos_var), line_align[0].clone());

                        // extract N surrounding bases around indel's position
                        // the resulted extracted sequence have N + INDEL_START + N
                        let n_surrounding: usize = 5;
                        let variant_flank = get_genomic_sequence(
                            chr_id,
                            cmp::max(indel.1 - n_surrounding, 0) as u64,
                            (indel.1 + n_surrounding) as u64,
                            &mut faidx_child,
                        );

                        let variant_char = indel.3;
                        // assert_eq!(pos_indels.len(), 1);

                        // Get hompolymer size around indel
                        let lh = expand_homopolymer(6, &[variant_char], &variant_flank);
                        // dbg!("HOMOPOLYMER SIZE", lh);

                        // Transform extracted sequence (flank) and indel first base (char) in &str
                        let variant_char_str = std::str::from_utf8(&[variant_char]).unwrap();
                        let variant_flank_str = std::str::from_utf8(&variant_flank).unwrap();

                        // Build VCF INFO object to store homopolymer length surrounding detected indel
                        // HS = Homopolymer length around potential incorrect indel;
                        // LTID = Locus_tag ids of interrupted CDS;
                        // EID = Entry ID that matched;
                        // IC = Indel bases context;
                        // BS = CDS start coordinate in initial genome;
                        // BE = CDS end coordinate in initial genome;
                        // ST = CDS strand;
                        let field = Field::new(
                            Key::Other("HS".parse().unwrap()),
                            Some(field::Value::Integer(lh as i32)),
                        );

                        let field_v2 = Field::new(
                            Key::Other("LTID".parse().unwrap()),
                            Some(field::Value::String(bf_names.join(",").parse().unwrap())),
                        );

                        let field_v3 = Field::new(
                            Key::Other("EID".parse().unwrap()),
                            Some(field::Value::String(first_bf_feature.sseqid.clone())),
                        );

                        let field_v4 = Field::new(
                            Key::Other("IC".parse().unwrap()),
                            Some(field::Value::String(variant_flank_str.to_string())),
                        );

                        let field_v5 = Field::new(
                            Key::Other("BS".parse().unwrap()),
                            Some(field::Value::Integer(genome_base_start as i32)),
                        );

                        let field_v6 = Field::new(
                            Key::Other("BE".parse().unwrap()),
                            Some(field::Value::Integer(genome_base_end as i32)),
                        );

                        let strand_char = if strand == 0 {
                            "+".to_string()
                        } else {
                            "-".to_string()
                        };

                        let field_v7 = Field::new(
                            Key::Other("ST".parse().unwrap()),
                            Some(field::Value::String(strand_char)),
                        );

                        let info = Info::try_from(vec![
                            field, field_v2, field_v3, field_v4, field_v5, field_v6, field_v7,
                        ])
                        .unwrap();

                        let ref_bases = match indel.0 {
                            // https://samtools.github.io/hts-specs/VCFv4.3.pdf - Examples in section 5
                            Ins => {
                                let base = get_genomic_sequence(
                                    chr_id,
                                    (indel.1 + 1) as u64,
                                    (indel.2) as u64,
                                    &mut faidx_child,
                                );

                                base
                            }
                            Del => {
                                let base = get_genomic_sequence(
                                    chr_id,
                                    (indel.1) as u64,
                                    (indel.2 - 1) as u64,
                                    &mut faidx_child,
                                );
                                let del_size = indel.2 - indel.1;
                                let duplicated_base: Vec<_> =
                                    base.into_iter().cycle().take(del_size + 1).collect();

                                duplicated_base
                                // let variant_char = if strand == 0 {
                                //     get_genomic_sequence(
                                //         chr_id,
                                //         (indel.1 + 1) as u64,
                                //         (indel.2 + 1) as u64,
                                //         &mut faidx_child,
                                //     )
                                // } else {
                                //     get_genomic_sequence(
                                //         chr_id,
                                //         (indel.1) as u64,
                                //         (indel.2) as u64,
                                //         &mut faidx_child,
                                //     )
                                // };
                                // variant_char
                            }
                            _ => panic!(
                                "Only variant of types insertion and deletion are considered"
                            ),
                        };

                        let alt_bases = match indel.0 {
                            Ins => {
                                let base = get_genomic_sequence(
                                    chr_id,
                                    (indel.1 + 1) as u64,
                                    (indel.2) as u64,
                                    &mut faidx_child,
                                );
                                let del_size = indel.2 - indel.1;
                                let duplicated_base: Vec<_> =
                                    base.into_iter().cycle().take(del_size + 1).collect();

                                duplicated_base
                                // let variant_char = if strand == 0 {
                                //     get_genomic_sequence(
                                //         chr_id,
                                //         (indel.1 + 1) as u64,
                                //         (indel.2 + 1) as u64,
                                //         &mut faidx_child,
                                //     )
                                // } else {
                                //     get_genomic_sequence(
                                //         chr_id,
                                //         (indel.1) as u64,
                                //         (indel.2) as u64,
                                //         &mut faidx_child,
                                //     )
                                // };
                                // variant_char
                            }
                            Del => [variant_char.clone()].to_vec(),

                            // get_genomic_sequence(
                            //     chr_id,
                            //     (indel.1 + 1) as u64,
                            //     (indel.2) as u64,
                            //     &mut faidx_child,
                            // ),
                            _ => panic!(
                                "Only variant of types insertion and deletion are considered"
                            ),
                        };

                        // Transform extracted sequence (flank) and indel first base (char) in &str
                        let ref_bases_str = std::str::from_utf8(&ref_bases).unwrap();
                        let alt_bases_str = std::str::from_utf8(&alt_bases).unwrap();

                        let record = vcf::Record::builder()
                            // .set_chromosome(chr_id.parse().unwrap())
                            .set_chromosome(chr_id.parse().unwrap())
                            .set_position(Position::from(pos_var))
                            .set_reference_bases(alt_bases_str.parse().unwrap())
                            .set_alternate_bases(ref_bases_str.parse().unwrap())
                            .set_info(info)
                            .build()
                            .unwrap();

                        let mut vcf_records_to_write =
                            vcf_records_to_write_child.lock().expect("Erro no vcf");
                        vcf_records_to_write.push(record);
                    }
                } else {
                    let mut entries_no_dna = entries_no_dna_child
                        .lock()
                        .expect("Error writing parallel entries names.");
                    entries_no_dna.push(query.to_string());
                }

                // .filter_map(|(mut index, s)| {
                //     // let mut altbases = String::new();
                //     // altbases.push(ena.seq()[alignment.ystart + index] as char);
                //     // altbases.push(boundary_seq[index] as char);
                //     if alignment.ystart + index >= ena.seq().len() {
                //         return None
                //     }

                //     if index >= boundary_seq.len() {
                //         return None
                //     }

                //     Some(
                //         (
                //         // we're dealing with glocal alignment, therefore, the sequecence y
                //         // can start after first base
                //         ena.id(),
                //         alignment.ystart + index,
                //         s,
                //         (boundary_seq[index] as char).to_string(),
                //         (ena.seq()[alignment.ystart + index] as char).to_string(),
                //         // boundary_seq_start + index as u64,
                //         (alignment.xstart + index) as u64 + boundary_seq_start - 1,
                //         alignment.xstart + index,
                //         )
                //     )
                // })

                // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

                // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {

                // let ena_record = reader.records().next().unwrap();
            };

            let mut output_fasta = output_fasta_child.lock().expect("Error on fasta");
            output_fasta.push(record);
        });

    // for (bf1, bf2) in count_possible_frameshifts.iter() {
    // if let Some(query) = bf1.sseqid.strip_prefix("UniRef90_") {
    //     let query_url = format!(
    //         "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
    //         // "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
    //         query
    //     );

    //     let req = ureq::get(&query_url)
    //         .set("Accept", "application/xml")
    //         .call();

    //     let reader = match req {
    //         Ok(r) => get_record_from_ena(r),
    //         Err(_e) => {
    //             let query_url = format!(
    //                 // "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
    //                 "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
    //                 query
    //             );

    //             let req = ureq::get(&query_url)
    //                 .set("Accept", "application/xml")
    //                 .call()
    //                 .unwrap();

    //             eprintln!(
    //                 "Error when mapping to reciprocal ProteinKB entry for UniRef90_{}",
    //                 &query
    //             );
    //             get_record_from_insilico_translation(req, query)
    //         }
    //     };
    //     // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

    //     // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {

    //     match hash_cds_fasta.get(&bf1.qseqid) {
    //         Some(record) => output_fasta.push(record.clone()),
    //         None => println!("Error, couldnt find sequence"),
    //     }

    //     match hash_cds_fasta.get(&bf2.qseqid) {
    //         Some(record) => output_fasta.push(record.clone()),
    //         None => println!("Error, couldnt find sequence"),
    //     }

    //     // let ena_record = reader.records().next().unwrap();

    //     let ena = match reader {
    //         Ok(rec) => {
    //             output_fasta.push(rec.clone());
    //             rec
    //         }
    //         Err(_e) => panic!("Could not find Record"),
    //     };

    //     let ena_seq_str = std::str::from_utf8(ena.seq()).unwrap();

    // if let Some(query) = bf1.sseqid.strip_prefix("UniRef90_") {
    //     let query_url = format!(
    //         "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
    //         // "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
    //         &query
    //     );

    //     println!("{}", &query_url);

    //     let visual_url = format!("https://www.uniprot.org/uniprot/{}", &query);

    //     debug!(
    //         "Performing an HTTP GET request against {} to retrieve DNA sequence for entry {}.",
    //         &visual_url, &bf1.sseqid
    //     );
    //     let req = ureq::get(&query_url)
    //         .set("Accept", "application/xml")
    //         .call();

    //     let reader = match req {
    //         Ok(r) => get_record_from_ena(r),
    //         Err(_e) => {
    //             let query_url = format!(
    //                 // "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
    //                 "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
    //                 query
    //             );

    //             let req = ureq::get(&query_url)
    //                 .set("Accept", "application/xml")
    //                 .call()
    //                 .unwrap();

    //             warn!(
    //                 "Could not retrieve DNA sequence from entry UniRef90_{}",
    //                 &query
    //             );
    //             // get_record_from_insilico_translation(req, query)
    //             unimplemented!()
    //         }
    //     };
    //     // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

    //     // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {

    //     // let ena_record = reader.records().next().unwrap();

    //     let ena_rec = match reader {
    //         Ok(rec) => {
    //             // output_fasta.push(rec.clone());
    //             Ok(rec)
    //         }
    //         Err(e) => {
    //             eprintln!("Could not find Record");
    //             Err(e)
    //         }
    //     };
    // };

    // let (tbl_bf1, chr_id) = hash_tbl_feat.get(&bf1.qseqid).unwrap();
    // let (tbl_bf2, chr_id2) = hash_tbl_feat.get(&bf2.qseqid).unwrap();

    // assert_eq!(chr_id, chr_id2);
    // let boundary_seq_description;
    // let boundary_seq = if tbl_bf1.start < tbl_bf2.start {
    //     boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf1.start, tbl_bf2.end);
    //     get_genomic_sequence(chr_id, tbl_bf1.start, tbl_bf2.end, &mut faidx)
    // } else {
    //     boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf2.start, tbl_bf1.end);
    //     get_genomic_sequence(chr_id, tbl_bf2.start, tbl_bf1.end, &mut faidx)
    // };

    // let read_id = format!("my_genomic_boundary_{}_{}", bf1.qseqid, bf2.qseqid);
    // let description = Some(&boundary_seq_description[..]);
    // // let sequence = b"ACGT";
    // let record = Record::with_attrs(&read_id[..], description, &boundary_seq);
    // // output_fasta.push(record);

    // let boundary_seq_str = std::str::from_utf8(&boundary_seq).unwrap();

    // let _chunks_diff = dissimilar::diff(ena_seq_str, boundary_seq_str);

    // println!("CHUNKS = {:?}", &chunks_diff);
    // ... process the Uniprot entry ...
    // }
    // } else {
    //     panic!("Error: protein subject is not from UniRef90 database.")
    // }
    // }

    // Uncomment the following snippet if you would like to stdout fasta genomic sequence of boundary seq among joined features
    // START SNIPPET
    // let mut writer = fasta::Writer::new(std::io::stdout());
    // let output_fasta_ref_to_write = output_fasta_ref.lock().unwrap();
    // for r in output_fasta_ref_to_write.iter() {
    //     writer
    //         .write_record(r)
    //         .expect("Error while writing FASTA output");
    // }
    // END SNIPPET

    let mut vcf_records_to_write = vcf_records_to_write_ref.lock().unwrap();

    vcf_records_to_write.sort_by_key(|r: &vcf::Record| (r.chromosome().to_string(), r.position()));

    let mut header = vcf::Header::builder();

    for vcf_rec in vcf_records_to_write.iter() {
        let chr_id = vcf_rec.chromosome().to_string();
        // if !contigs.contains(&chr_id) {
        let contig = Map::<Contig>::new(chr_id.parse().unwrap());
        header = header.add_contig(contig.clone());
        // .build();

        // contigs.insert(chr_id);
        // }
    }

    // let meta = header::record::value::Map::<header::record::value::map::Meta>::new(
    //     "Assay",
    //     vec![String::from("WholeGenome"), String::from("Exome")],
    // );

    let map_field = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("HS".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::Integer,
        "Homopolymer length around potential incorrect indel",
    );

    let map_field_v2 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("LTID".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::String,
        "Locus_tag ids of interrupted CDS",
    );

    let map_field_v3 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("EID".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::String,
        "Entry ID that matched",
    );

    let map_field_v4 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("IC".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::String,
        "Indel bases context",
    );

    let map_field_v5 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("BS".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::Integer,
        "CDS start coordinate in initial genome",
    );

    let map_field_v6 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("BE".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::Integer,
        "CDS end coordinate in initial genome",
    );

    let map_field_v7 = header::record::value::Map::<header::record::value::map::Info>::new(
        Key::Other("ST".parse().unwrap()),
        header::Number::Count(1),
        header::info::Type::String,
        "CDS strand",
    );

    let header_final = header
        .add_info(map_field)
        .add_info(map_field_v2)
        .add_info(map_field_v3)
        .add_info(map_field_v4)
        .add_info(map_field_v5)
        .add_info(map_field_v6)
        .add_info(map_field_v7)
        .build();

    writer_vcf
        .write_header(&header_final)
        .expect("Couldn't write VCF file contig header.");

    let hash_alignment_pretty_to_write = hash_alignment_pretty_ref.lock().unwrap();
    for vcf_rec in vcf_records_to_write.iter() {
        let var_position = usize::from(vcf_rec.position());
        let var_ref = vcf_rec.reference_bases();
        let var_alt = vcf_rec.alternate_bases();
        let chr_key = vcf_rec.chromosome().to_string();
        let lines_aln = hash_alignment_pretty_to_write
            .get(&(chr_key.clone(), var_position))
            .unwrap();
        info!("===============================================");
        info!(
            "Var genomic pos: {}:{} | REF: {} | ALT: {}",
            &chr_key, var_position, var_ref, var_alt
        );
        for l in lines_aln.split('\n') {
            info!("{}", l);
        }
        info!("===============================================");

        writer_vcf
            .write_record(vcf_rec)
            .expect("Couldnt write record.");
    }

    info!("OVERALL STATS");
    info!("===============================================");
    info!("Num input sequences ....................................:");
    info!(
        "Num consecutive features (CDS) matching same best hit (possible frameshift) .: {}",
        &count_possible_frameshifts_num
    );
    info!(
        "Num grouped consecutive features (CDS) associated with frameshifts .: {}",
        &count_possible_frameshifts_v2_num
    );
    let pair_with_indel_all = pair_with_indel_ref.lock().unwrap();
    info!(
        "Num feature pairs associated with INDELs ......: {}",
        &pair_with_indel_all
    );
    info!(
        "Num homopolymers associated with INDELs ......: {}",
        &count_possible_frameshifts_v2_num
    );
    let entries_no_dna_all = entries_no_dna_ref.lock().unwrap();
    info!(
        "Num entries which doesnt have genomic dna associated .: {}",
        &entries_no_dna_all.len()
    );
    info!(
        "Num complex loci (not a homopolymer) ..................................: {}",
        &complex_loci_ref.lock().unwrap()
    );
    info!("Num actions ..................................:");
    info!("Num insertions ...............................:");
    info!("Num deletions ................................:");

    info!("Results in out.vcf.gz",);

    // let query_url = format!(
    //     "https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+{}&format=xml&compress=yes",
    //     query
    // );

    // let req = ureq::get(&query_url).set("Accept", "application/xml");
    // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

    // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {
    //     let entry = r.unwrap();
    //     for db_ref in entry.db_references.iter() {
    //         for prop in db_ref.property.iter() {
    //             if prop.value == "Genomic_DNA" {
    //                 dbg!(&db_ref);
    //             }
    //         }
    //     }
    //     // ... process the Uniprot entry ...
    // }

    // // let query_url_ena = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases?style=json";
    // let query_url_ena =
    //     "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/ena_sequence/L12344,L12345/fasta";
    // let req_ena = ureq::get(&query_url_ena);

    // let reader_ena = fasta::Reader::new(req_ena.call().unwrap().into_reader());

    // for result in reader_ena.records() {
    //     dbg!(result);
    // }

    // dbg!(reader_ena);
}
