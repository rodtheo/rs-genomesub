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
use flate2::write::GzEncoder;
use flate2::Compression;
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
use std::collections::HashMap;

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
    self as vcf,
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
use std::sync::{Arc, Mutex};

/// Given a Protein ID, get the corresponding genomic sequence through an API request to ENA server.
pub fn get_record_from_ena(response: Response) -> Result<Record, std::io::Error> {
    let reader = libflate::gzip::Decoder::new(response.into_reader()).unwrap();

    let entry = uniprot::uniprot::parse(std::io::BufReader::new(reader))
        .next()
        .unwrap();

    let mut prop_ena_id: Option<String> = None;
    for db_ref in entry.unwrap().db_references.iter() {
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

    reader_ena.records().next().unwrap()
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
    while i_threeprime < reference.len() && ([reference[pos + i_threeprime]] == *alt) {
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

    let input_path_tbl_str = matches.value_of("input").unwrap();
    let input_path_genome_str = matches.value_of("genome").unwrap();
    let input_path_genes_str = matches.value_of("genes").unwrap();
    let input_path_diamond_str = matches.value_of("diamond").unwrap();

    let input_path = PathBuf::from(&input_path_genes_str);

    let mut input = File::open(input_path).unwrap();
    let res = fasta::Reader::new(&mut input);

    let out_vcf = File::create("out.vcf.gz")
        .map(|x| GzEncoder::new(x, Compression::default()))
        .expect("Unable to create file to write output sequences.");

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

    // let mut output_fasta: Vec<fasta::Record> = Vec::new();
    let mut output_fasta_ref = Arc::new(Mutex::new(Vec::new()));
    let mut vcf_records_to_write_ref = Arc::new(Mutex::new(Vec::new()));

    let input_path_diamond = PathBuf::from(&input_path_diamond_str);
    let fblastout = fs::read_to_string(input_path_diamond).expect("Unable to read diamond result");
    // let fblastout = include_str!("../../examples/output_diamond_sample.txt");

    // Parse diamond result into BlastFeature struct
    let (_inp, res) = parse_input(&fblastout).unwrap();

    let input_path_genome = PathBuf::from(&input_path_genome_str);
    let path_fsa = Path::new(&input_path_genome);

    let mut faidx = IndexedReader::from_file(&path_fsa).unwrap();

    std::env::set_var("RAYON_NUM_THREADS", "4");

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

    eprintln!(
        "Quantity of pairs of annotated coding sequences (CDS) matching the same uniref entry: {}.\nThis fact suggests that there may be frameshifts interrupting a pair.",
        &count_possible_frameshifts.len()
    );

    eprintln!(
        "V2: Quantity of pairs of annotated coding sequences (CDS) matching the same uniref entry: {}.\nThis fact suggests that there may be frameshifts interrupting a pair.",
        &count_possible_frameshifts_v2.len()
    );

    // let mut boundary_seq_description_v2 = format!("{}:{}-{}", chr_id, tbl_bf1.start, tbl_bf2.end);

    count_possible_frameshifts_v2
        .par_iter_mut()
        .for_each(|pfs| {
            let mut faidx_child = IndexedReader::from_file(&path_fsa).unwrap();
            let output_fasta_child = output_fasta_ref.clone();
            let vcf_records_to_write_child = vcf_records_to_write_ref.clone();

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
                    Err(_e) => {
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
                        unimplemented!()
                    }
                };

                let ena_rec = match reader {
                    Ok(rec) => {
                        // output_fasta.push(rec.clone());
                        Ok(rec)
                    }
                    Err(e) => {
                        eprintln!("Could not find Record");
                        Err(e)
                    }
                };

                let ena = ena_rec.unwrap();

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
                        let var_pos = idx + genome_base_start as usize;
                        (op, var_pos, var_pos + 1)
                    })
                    .collect::<Vec<(&AlignmentOperation, usize, usize)>>();

                // Iterate over each indel observed in best alignment
                for indel in pos_indels.iter() {
                    // extract N surrounding bases around indel's position
                    // the resulted extracted sequence have N + INDEL_START + N
                    let n_surrounding: usize = 5;
                    let variant_flank = get_genomic_sequence(
                        chr_id,
                        cmp::max(indel.1 - n_surrounding, 0) as u64,
                        (indel.1 + n_surrounding) as u64,
                        &mut faidx_child,
                    );

                    // assert_eq!(pos_indels.len(), 1);

                    let strand = best_aln.0;
                    // extract indel first base
                    // strand == 0 -> Forward
                    let variant_char = if strand == 0 {
                        get_genomic_sequence(
                            chr_id,
                            (indel.1 + 1) as u64,
                            (indel.2) as u64,
                            &mut faidx_child,
                        )
                    } else {
                        get_genomic_sequence(
                            chr_id,
                            (indel.1) as u64,
                            (indel.2 - 1) as u64,
                            &mut faidx_child,
                        )
                    };

                    dbg!(&variant_char);

                    // Get hompolymer size around indel
                    let lh = expand_homopolymer(6, &variant_char, &variant_flank);
                    dbg!("HOMOPOLYMER SIZE", lh);

                    // Transform extracted sequence (flank) and indel first base (char) in &str
                    let variant_char_str = std::str::from_utf8(&variant_char).unwrap();
                    let variant_flank_str = std::str::from_utf8(&variant_flank).unwrap();

                    // Indel start position
                    let pos_var = (indel.1 + 1) as usize;

                    // Build VCF INFO object to store homopolymer length surrounding detected indel
                    let field = Field::new(
                        Key::SamplesWithDataCount,
                        Some(field::Value::Integer(lh as i32)),
                    );

                    let field_v2 = Field::new(
                        Key::Other("PP".parse().unwrap()),
                        Some(field::Value::String(bf_names.join("_").parse().unwrap())),
                    );

                    let field_v3 = Field::new(
                        Key::Other("UR".parse().unwrap()),
                        Some(field::Value::String(first_bf_feature.sseqid.clone())),
                    );

                    let field_v4 = Field::new(
                        Key::Other("AF".parse().unwrap()),
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
                                (indel.1 + 1) as u64,
                                (indel.2) as u64,
                                &mut faidx_child,
                            );
                            let duplicated_base: Vec<_> =
                                base.into_iter().cycle().take(2).collect();

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
                        _ => panic!("Only variant of types insertion and deletion are considered"),
                    };

                    let alt_bases = match indel.0 {
                        Ins => {
                            let base = get_genomic_sequence(
                                chr_id,
                                (indel.1 + 1) as u64,
                                (indel.2) as u64,
                                &mut faidx_child,
                            );
                            let duplicated_base: Vec<_> =
                                base.into_iter().cycle().take(2).collect();

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
                        Del => get_genomic_sequence(
                            chr_id,
                            (indel.1 + 1) as u64,
                            (indel.2) as u64,
                            &mut faidx_child,
                        ),
                        _ => panic!("Only variant of types insertion and deletion are considered"),
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

                    if first_bf_feature.sseqid == "UniRef90_A0A0H3KZC3" {
                        // dbg!(best_aln);
                        dbg!(best_aln.1.operations.len());
                        // dbg!(pos_indels);
                    }
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

    let mut writer = fasta::Writer::new(std::io::stdout());
    let output_fasta_ref_to_write = output_fasta_ref.lock().unwrap();
    for r in output_fasta_ref_to_write.iter() {
        writer
            .write_record(r)
            .expect("Error while writing FASTA output");
    }

    let mut vcf_records_to_write = vcf_records_to_write_ref.lock().unwrap();

    vcf_records_to_write.sort_by_key(|r: &vcf::Record| r.position());

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

    let header_final = header.build();
    writer_vcf
        .write_header(&header_final)
        .expect("Coldnt write VCF file contig header.");

    for vcf_rec in vcf_records_to_write.iter() {
        writer_vcf.write_record(&vcf_rec);
    }

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
