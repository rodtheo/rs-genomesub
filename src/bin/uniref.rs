extern crate libflate;
extern crate uniprot;
extern crate ureq;

use multipeek::multipeek;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;
use bio::alphabets::dna;
use bio::io::fasta::{self, IndexedReader, Record};
use bio::alignment::sparse::*;
use noodles::vcf::{self as vcf, header::Contig, record::reference_bases::Base, record::Position};
// use nom::error::dbg_dmp;
// use rs_genomesub::tbl::Feature;
// use uniprot::uniprot::DbReference;
use ureq::Response;

extern crate rs_genomesub;
use crate::rs_genomesub::blast_tabular::*;
use crate::rs_genomesub::tbl;

use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::HashMap;

use std::fs;
use std::fs::File;
// use std::io::Write;
use std::path::Path;
// use std::path::PathBuf;

use clap::{arg, Arg, Command};
use min_max::*;
use rayon::prelude::*;
use simplelog::*;
use std::collections::HashSet;
use std::error::Error;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use protein_translate::translate;

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct BoundarySeq {
    contig: String,
    start_assembly: u64,
    end_assembly: u64,
    boundary_cor_seq: Vec<u8>,
}

/// Given a Protein ID, get the corresponding genomic sequence through an API request to ENA server.
pub fn get_record_from_ena(response: Response) -> Result<Record, Box<dyn Error>> {
    let reader = libflate::gzip::Decoder::new(response.into_reader()).unwrap();

    let entry = uniprot::uniprot::parse(std::io::BufReader::new(reader))
        .next()
        .unwrap();

    let r = match entry {
        Ok(r) => r,
        Err(e) => return Err(Box::new(e)),
    };

    let mut prop_ena_id: Option<String> = None;
    for db_ref in r.db_references.iter() {
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

    match reader_ena.records().next().unwrap() {
        Ok(re) => return Ok(re),
        Err(e) => return Err(Box::new(e)),
    }
    // reader_ena.records().next().unwrap()
}

pub fn get_record_from_insilico_translation(
    response: Response,
    prot_id: &str,
) -> Result<Record, Box<dyn Error>> {
    let reader = libflate::gzip::Decoder::new(response.into_reader()).unwrap();

    let entry = uniprot::uniref::parse(std::io::BufReader::new(reader))
        .next()
        .unwrap();

    let protseq = if let uniprot::uniref::Member { sequence: seq, .. } =
        entry.unwrap().representative_member
    {
        if let Some(uniprot::uniparc::Sequence {
            sequence: string_sequence,
            ..
        }) = seq
        {
            string_sequence
        } else {
            unimplemented!()
        }
    } else {
        unimplemented!()
    };

    // eprintln!("STRING SEQ: {}", &protseq);

    // println!("{:?}", &entry);

    // let protseq = b"ATC";
    // let read_id = "ehissoai";
    let read_id = format!("InSilico_{}", prot_id);
    let protseq_u8 = protseq.as_bytes();
    let dna = translate(protseq_u8);

    let description = None;
    let record = Record::with_attrs(&read_id, description, dna.as_bytes());

    Ok(record)
}

pub fn apply_variant(
    boundary_seq: & Vec<u8>,
    variant_pos: usize,
    alternate_bases: String,
    mut_type: &AlignmentOperation,
    is_revcomp: bool,
) -> Vec<u8> {
    match mut_type {
        Ins => {
            let seq = [
                &boundary_seq[..variant_pos],
                &boundary_seq[variant_pos + alternate_bases.len() - 1..],
            ]
            .concat();
            if is_revcomp {
                dna::revcomp(seq)
            } else {
                seq
            }
        }

        Del => {
            let seq = [
                &boundary_seq[..variant_pos],
                &alternate_bases.as_bytes().to_vec(),
                &boundary_seq[variant_pos + 1..],
            ]
            .concat();

            if is_revcomp {
                dna::revcomp(seq)
            } else {
                seq
            }
        }
        _ => Vec::<u8>::new(),
    }
}

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
        .about("Inspect candidate frame shifts errors within coding regions")
        .arg_required_else_help(true)
        .arg(
            Arg::new("tbl_in")
                .value_name("TBL_IN")
                .help("Specify the annotations in tbl format.")
                .required(true), // .min_values(1),
        )
        .arg(arg!(-a --assembly <ASSEMBLY_IN>).required(true).help("Assembly in FASTA (requires faidx index)"))
        // .arg(
        //     Arg::new("genome")
        //         .value_name("GENOME")
        //         .short('r')
        //         .help("Assembly in FASTA. Requires faidx index.")
        //         .takes_value(true),
        // )
        .arg(
            Arg::new("genes_in")
                .value_name("GENES_IN")
                .help("Specify the file containing gene sequences (nucleotide FASTA). Make sure that sequence id matches the locus_tag in feature table (tbl) input.")
                .takes_value(true)
                .short('g')
                .long("genes")
                .required(true),
        )
        .arg(
            Arg::new("diamond_in")
                .value_name("DIAMOND_IN")
                .help("Specify the input diamond result obtained after running blastx with gene sequences")
                .takes_value(true)
                .short('b')
                .long("diamond")
                .required(true),
        )
        .arg(
            Arg::new("assembly_out")
                .value_name("ASSEMBLY_OUT")
                .help("Specify the output FASTA file name.")
                .takes_value(true)
                .short('o')
                .long("out-assembly")
                .required(true),
        )
        .arg(arg!(
            -t --threads ... "Set the number of threads."
        ).takes_value(true).default_value("4"))
        .arg(arg!(
            -d --debug ... "Set verbosity to debug level."
        ).takes_value(false))
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


    let input_path_tbl_str = matches.value_of("tbl_in").unwrap();
    let input_path_genome_str = matches.value_of("assembly").unwrap();
    let input_path_genes_str = matches.value_of("genes_in").unwrap();
    let input_path_diamond_str = matches.value_of("diamond_in").unwrap();
    let output_path_assembly_str = matches.value_of("assembly_out").unwrap();
    let n_threads = matches.value_of("threads").unwrap();
    
    let output_path_genome = std::path::Path::new(&output_path_assembly_str);
    let prefix = output_path_genome.parent().unwrap();
    std::fs::create_dir_all(prefix).unwrap();

    let mut out_log = PathBuf::from(&output_path_assembly_str);
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

    

    std::env::set_var("RAYON_NUM_THREADS", n_threads);

    let input_path = PathBuf::from(&input_path_genes_str);

    let mut input = File::open(input_path).unwrap();
    let res = fasta::Reader::new(&mut input);

    let mut hash_cds_fasta = HashMap::new();

    let mut record;
    let mut record_id;
    for s in res.records() {
        record = s.unwrap();
        record_id = record.id().to_string();
        hash_cds_fasta.insert(record_id, record);
    }

    info!(
        "Reading {} genes sequences from {}.",
        hash_cds_fasta.len(),
        &input_path_genes_str
    );

    let input_path_tbl = PathBuf::from(&input_path_tbl_str);
    let mut input_tbl = File::open(input_path_tbl).unwrap();
    let res_tbl = tbl::Reader::new(&mut input_tbl);
    let mut hash_tbl_feat = HashMap::new();

    for tbl_entry in res_tbl.iter() {
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

    let input_path_diamond = PathBuf::from(&input_path_diamond_str);
    let fblastout = fs::read_to_string(input_path_diamond).expect("Unable to read diamond result");
    // let fblastout = include_str!("../../examples/output_diamond_sample.txt");

    // Parse diamond result into BlastFeature struct
    let (_inp, res) = parse_input(&fblastout).unwrap();

    let input_path_genome = PathBuf::from(&input_path_genome_str);
    // let mut output_path_genome = PathBuf::from(&input_path_genome_str);
    let path_fsa = Path::new(&input_path_genome);

    // let mut faidx = IndexedReader::from_file(&path_fsa).unwrap();
    // let mut faidx = Arc::new(Mutex::new(IndexedReader::from_file(&path_fsa).unwrap()));
    // let faidx = Arc::new(IndexedReader::from_file(&path_fsa).unwrap());
    // warn!("OUTPUT PATH = {}",&output_path_assembly_str);
    
    
    // output_path_genome.set_file_name("out_cor.fasta");
    let output_assembly = File::create(output_path_genome).unwrap();
    // let ofile = File::create(output_path_genome).expect("Unable to create file to write output sequences.");
    // let mut writer_vcf = vcf::Writer::new(ofile);

    // let count_possible_frameshifts = res
    //     .into_iter()
    //     .tuple_windows()
    //     .filter(|(s1, s2)| s1.sseqid == s2.sseqid)
    //     .count();

    info!("Starting <green>{} workers</>.", n_threads);
    info!("Detecting adjacent genes on the same strand that give hits against the same subject (common BLAST hit).");
    let count_possible_frameshifts: Vec<(BlastFeature, BlastFeature)> = res
        // .into_iter()
        // .tuple_windows()
        .par_windows(2)
        .filter_map(|s1| {
            if (s1[0].sseqid == s1[1].sseqid) {
                Some((s1[0].clone(), s1[1].clone()))
            } else {
                None
            }
        })
        // .filter(|[s1, s2]| (s1.sseqid == s2.sseqid))
        // .map(|[s1, s2]| (s1, s2))
        .collect();

    let cfsorted = count_possible_frameshifts.clone();
    cfsorted.clone().sort_by(|a, b| {
        let a_start = hash_tbl_feat[&a.0.qseqid].0.start;
        let b_start = hash_tbl_feat[&b.0.qseqid].0.start;
        if (a_start < b_start) {
            Ordering::Less
        } else if (a_start == b_start) {
            Ordering::Equal
        } else {
            Ordering::Greater
        }
    });

    let mut grouped_frames = Vec::new();
    let mut iter_frames = multipeek(cfsorted.into_iter());    
    // let mut el_peek = iter_frames.peek().cloned();
    let mut el_peek = iter_frames.next();
    let mut inner_grouped = Vec::new();
    while !el_peek.is_none() {
        inner_grouped = Vec::new();
        let second_peek = iter_frames.peek_nth(0).cloned();
        let third_peek = iter_frames.peek_nth(1).cloned();

        let fel = el_peek.unwrap_or_default();
        let sel = second_peek.unwrap_or_default();
        let tel = third_peek.unwrap_or_default();

        // dbg!("BEGIN ==================");
        // dbg!(&fel);
        // dbg!(&sel);
        // dbg!(&tel);
        // dbg!("END ====================");
        

        // if fel.0.sseqid == fel.1.sseqid {
        //     inner_grouped.push(fel.clone())
        // }

        if fel.0.sseqid == sel.0.sseqid && sel.0.sseqid != tel.0.sseqid {
            
            inner_grouped.push(fel);
            inner_grouped.push(sel);
            
            iter_frames.next();
            el_peek = iter_frames.next();
        }
        else if fel.0.sseqid == sel.0.sseqid && sel.0.sseqid == tel.0.sseqid {
         
            inner_grouped.push(fel);
            inner_grouped.push(sel);
            inner_grouped.push(tel);

            iter_frames.next();
            iter_frames.next();
            el_peek = iter_frames.next();
        } else {
            inner_grouped.push(fel);
            
            el_peek = iter_frames.next();
        }

        // if inner_grouped.len() > 0 {
            grouped_frames.push(inner_grouped);
        // }
        
    }
    
    // let fg = grouped_frames.iter().map(|x| x)
                // .filter(|v| v.len() > 1);
    // dbg!(&grouped_frames);
    // dbg!(count_pos);

    // let count_possible_frameshifts = res
    //     .into_iter()
    //     .tuple_windows()
    //     .filter(|(s1, s2)| (s1.sseqid == s2.sseqid))
    //     .collect::<Vec<(BlastFeature, BlastFeature)>>();

    info!(
        "{} Grouped-pairs of adjacent genes matching same blast hit",
        &grouped_frames.len()
    );

    info!(
        "{} pairs of adjacent genes matching same blast hit.",
        &count_possible_frameshifts.len()
    );

    // debug!("I can write <b>bold</b> text or use tags to <red>color it</>");

    let mut contigs = HashSet::new();
    let mut homology_loci = Arc::new(Mutex::new(HashSet::new()));
    let mut output_fasta_hash = Arc::new(Mutex::new(HashSet::new()));
    // let mut vcf_records_to_write: Vec<vcf::Record> = Vec::new();
    let mut vcf_records_to_write_ref = Arc::new(Mutex::new(Vec::new()));
    let mut seqs_to_change_ref = Arc::new(Mutex::new(Vec::new()));

    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(4)
    //     .build_global()
    //     .unwrap();

    grouped_frames
        .par_iter()
        .enumerate()
        .for_each(|(i, blast_features)| {
            // https://stackoverflow.com/questions/30559073/cannot-borrow-captured-outer-variable-in-an-fn-closure-as-mutable
            // for (i, (bf1, bf2)) in count_possible_frameshifts.iter().enumerate() {
            // let faidx_child = faidx.clone();
            let output_fasta_child = output_fasta_ref.clone();
            let vcf_records_to_write_child = vcf_records_to_write_ref.clone();
            let seqs_to_change_child = seqs_to_change_ref.clone();
            let output_fasta_hash_child = output_fasta_hash.clone();
            let mut faidx_child = IndexedReader::from_file(&path_fsa).unwrap();
            let (bf1, bf2) = &blast_features[0];
            
            debug!(
                "Processing {} adjacent gene pair {}/{}: {} + {} that have matched {}.",
                &blast_features.len(),
                &i,
                &grouped_frames.len(),
                &bf1.qseqid,
                &bf2.qseqid,
                &bf1.sseqid
            );

            if let Some(query) = bf1.sseqid.strip_prefix("UniRef90_") {
                let query_url = format!(
                    "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
                    // "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
                    &query
                );

                let visual_url = format!("https://www.uniprot.org/uniprot/{}", &query);

                debug!(
                    "Performing an HTTP GET request against {} to retrieve DNA sequence for entry {}.",
                    &visual_url, &bf1.sseqid
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
                        get_record_from_insilico_translation(req, query)
                    }
                };
                // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

                // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {

                // let ena_record = reader.records().next().unwrap();

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

                if let Ok(ena) = ena_rec {
                    // let ena_seq_str = std::str::from_utf8(ena.seq()).unwrap();

                    let (tbl_bf1, chr_id) = hash_tbl_feat
                        .get(&bf1.qseqid)
                        .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));
                    let (tbl_bf2, chr_id2) = hash_tbl_feat
                        .get(&bf2.qseqid)
                        .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));

                    assert_eq!(chr_id, chr_id2);

                    

                    // if !contigs.contains(chr_id) {
                    //     let header = vcf::Header::builder()
                    //         .add_contig(Contig::new(chr_id))
                    //         .build();
                    //     writer_vcf
                    //         .write_header(&header)
                    //         .expect("Coldnt write VCF file contig header.");
                    //     contigs.insert(chr_id);
                    // }


                    let mut boundary_seq_start =
                        blast_features.clone().into_iter().fold(Vec::new(), |mut acc, (bf1, bf2)| {
                            let (tbl_bf1, chr_id) = hash_tbl_feat
                                .get(&bf1.qseqid)
                                .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));
                            let (tbl_bf2, chr_id2) = hash_tbl_feat
                                .get(&bf2.qseqid)
                                .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));
                            acc.push(tbl_bf1.start);
                            acc.push(tbl_bf2.start);
                            acc}).into_iter().min().unwrap();
                    
                        // min!(tbl_bf1.start, tbl_bf2.start, tbl_bf1.end, tbl_bf2.end);
                    let mut boundary_seq_end = blast_features.iter().fold(Vec::new(), |mut acc_n, (bf1, bf2)| {
                        let (tbl_bf1, chr_id) = hash_tbl_feat
                            .get(&bf1.qseqid)
                            .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));
                        let (tbl_bf2, chr_id2) = hash_tbl_feat
                            .get(&bf2.qseqid)
                            .expect(&format!("Gene {} is note a Key in gene dict", &bf1.qseqid));
                        acc_n.push(tbl_bf1.end);
                        acc_n.push(tbl_bf2.end);
                        acc_n}).into_iter().max().unwrap();
                        
                    let boundary_seq_description =
                        format!("{}:{}-{}", chr_id, boundary_seq_start, boundary_seq_end);

                    // if blast_features.len() > 1 {
                    //         dbg!("########## VEC > 1");
                    //         dbg!(&boundary_seq_description);
                    //         dbg!(&blast_features);
                    //     }

                    // // Sparse alignment from kmer matches
                    // let mut chr_seq = Vec::new();
                    // faidx_child.fetch_all(chr_id)
                    //     .expect("Couldn't fetch interval");
                    // faidx_child.read(&mut chr_seq).expect("Coldn't read the interval");
                    
                    
                    // let qlen = ena.seq().len() as f64;
                    // let rlen = chr_seq.len() as f64;
                    // let klen = ((qlen + rlen).log2() / 2f64 ).ceil() as usize;

                    // let matches = find_kmer_matches(ena.seq(), &chr_seq, klen);
                    // let sparse_al = lcskpp(&matches, klen);
                    // let match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| matches[*i]).collect();
                    // // dbg!(match_path);
                    // let (start_x, start_y) = match_path.first().unwrap();
                    // let (end_x, end_y) = match_path.last().unwrap();
                    // dbg!(sparse_al.score);

                    let mut boundary_seq = get_genomic_sequence(
                        chr_id,
                        boundary_seq_start.clone(),
                        boundary_seq_end.clone(),
                        // start_y.clone() as u64,
                        // end_y.clone() as u64,
                        &mut faidx_child,
                    );

                    // let boundary_seq = if tbl_bf1.start < tbl_bf2.start {
                    //     boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf1.start, tbl_bf2.end);
                    //     get_genomic_sequence(chr_id, tbl_bf1.start, tbl_bf2.end, &mut faidx)
                    // } else {
                    //     boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf2.start, tbl_bf1.end);
                    //     get_genomic_sequence(chr_id, tbl_bf2.start, tbl_bf1.end, &mut faidx)
                    // };

                    // let read_id = format!("gb_{}_{}", bf1.qseqid, bf2.qseqid);
                    let ena_id = ena.id().split("|").last().unwrap();
                    let read_id = format!("gb_{}_{}", ena_id, i);
                    let description = Some(&boundary_seq_description[..]);
                    // let description = None;
                    // let sequence = b"ACGT";
                    let record_bondary =
                        Record::with_attrs(&read_id[..], description, &boundary_seq);

                    

                    // let boundary_seq_str = std::str::from_utf8(&boundary_seq).unwrap();

                    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
                    // gap open score: -5, gap extension score: -1

                    // x is the query or read sequence and y is the reference or template sequence.
                    let mut aligner =
                        Aligner::with_capacity(boundary_seq.len(), ena.seq().len(), -5, -1, &score);

                    let mut alignment = aligner.semiglobal(&boundary_seq, ena.seq());

                    let mut score_normalized =
                        
                        f64::from(alignment.score) / f64::from(boundary_seq.len() as i32);

                    let (match_count, mismatch_count) =
                        alignment.operations
                            .iter()
                            .fold((0, 0), |(matches, mismatches), x| match x {
                                Match => (matches + 1, mismatches),
                                Subst => (matches, mismatches + 1),
                                Del => (matches, mismatches),
                                Ins => (matches, mismatches),
                                Xclip(_) => (matches, mismatches),
                                Yclip(_) => (matches, mismatches),
                            });

                    let gap_free_for = (match_count as f64) / (match_count as f64 + mismatch_count as f64);

                    let boundary_seq_rev: Vec<u8> = dna::revcomp(&boundary_seq);

                    let alignment_rev = aligner.semiglobal(&boundary_seq_rev, ena.seq());

                    let score_normalized_rev =
                        f64::from(alignment_rev.score) / f64::from(boundary_seq_rev.len() as i32);
                    
                    let (match_count_rev, mismatch_count_rev) =
                        alignment_rev.operations
                            .iter()
                            .fold((0, 0), |(matches, mismatches), x| match x {
                                Match => (matches + 1, mismatches),
                                Subst => (matches, mismatches + 1),
                                Del => (matches, mismatches),
                                Ins => (matches, mismatches),
                                Xclip(_) => (matches, mismatches),
                                Yclip(_) => (matches, mismatches),
                            });
                    
                    let gap_free_rev = (match_count_rev as f64) / (match_count_rev as f64 + mismatch_count_rev as f64);
                    
                    warn!(
                        "Score normalized by query length FORWARD = {:.2}, REV = {:.2}",
                        &score_normalized, &score_normalized_rev,
                    );

                    warn!(
                        "GAP FREE Score normalized by query length FORWARD = {:.2}, REV = {:.2}",
                        &gap_free_for, &gap_free_rev,
                    );

                    let mut is_revcomp = false;

                    // if (score_normalized_rev >= 0.95 && score_normalized <= 0.95) {
                    if (gap_free_rev >= 0.95 && gap_free_for <= 0.95) {
                        alignment = alignment_rev;
                        boundary_seq = boundary_seq_rev;
                        is_revcomp = true;
                    }

                    // if (alignment.ylen > alignment.xlen + 30) {
                    //     dbg!("BEEEFORE ADJUSTING");
                    //     dbg!(&alignment.ylen);
                    //     dbg!(&alignment.xlen);
                    //     dbg!(&alignment.xstart);
                    //     dbg!(&alignment.ystart);
                    //     dbg!(&alignment.xend);
                    //     dbg!(&alignment.yend);
                    //     dbg!(&boundary_seq_start);
                    //     dbg!(&boundary_seq_end);

                    //     if is_revcomp {
                    //         dbg!("IS REV");
                    //         boundary_seq_end = boundary_seq_end.clone() + ((alignment.ystart - alignment.xstart) as u64);
                    //         boundary_seq_start = boundary_seq_start.clone() - ((alignment.ylen - alignment.xlen) as u64);
                    //     } else {
                    //         boundary_seq_end = boundary_seq_end.clone() + ((alignment.ylen - alignment.xlen) as u64);
                    //         boundary_seq_start = boundary_seq_start.clone() - ((alignment.ystart - alignment.xstart) as u64);
                    //     }


                    //     boundary_seq = get_genomic_sequence(
                    //         chr_id,
                    //         boundary_seq_start.clone(),
                    //         boundary_seq_end.clone(),
                    //         // start_y.clone() as u64,
                    //         // end_y.clone() as u64,
                    //         &mut faidx_child,
                    //     );

                    //     if is_revcomp {
                    //         boundary_seq = dna::revcomp(&boundary_seq);
                    //     }

                    //     // x is the query or read sequence and y is the reference or template sequence.
                    //     aligner =
                    //         Aligner::with_capacity(boundary_seq.len(), ena.seq().len(), -5, -1, &score);

                    //     alignment = aligner.semiglobal(&boundary_seq, ena.seq());

                    //     score_normalized =
                    //         f64::from(alignment.score) / f64::from(boundary_seq.len() as i32);


                    //         dbg!("AAAFTER ADJUSTING");
                    //         dbg!(&alignment.xstart);
                    //         dbg!(&alignment.ystart);
                    //         dbg!(&alignment.xend);
                    //         dbg!(&alignment.yend);
                    //         dbg!(&boundary_seq_start);
                    //         dbg!(&boundary_seq_end);
                    //         dbg!("Genes: {} + {} => {}", &bf1.qseqid, &bf2.qseqid, &ena.id());
                            
                    // }

                    // if (score_normalized >= 0.95 || score_normalized_rev >= 0.95) {
                    if (gap_free_for >= 0.95 || gap_free_rev >= 0.95) {
                        debug!("Genes: {} + {} => {}", &bf1.qseqid, &bf2.qseqid, &ena.id());
                        // debug!("{}", &alignment.score);
                        // debug!("{}", boundary_seq.len());
                        
                        let pos = alignment
                            .operations
                            .iter()
                            .enumerate()
                            .filter(|(_, &s)| match s {
                                Ins => true,
                                Del => true,
                                _ => false,
                            })
                            .filter_map(|(mut index, s)| {
                                // let mut altbases = String::new();
                                // altbases.push(ena.seq()[alignment.ystart + index] as char);
                                // altbases.push(boundary_seq[index] as char);
                                if alignment.ystart + index >= ena.seq().len() {
                                    return None
                                }

                                if index >= boundary_seq.len() {
                                    return None
                                }

                                Some(
                                    (
                                    // we're dealing with glocal alignment, therefore, the sequecence y
                                    // can start after first base
                                    ena.id(),
                                    alignment.ystart + index,
                                    s,
                                    (boundary_seq[index] as char).to_string(),
                                    (ena.seq()[alignment.ystart + index] as char).to_string(),
                                    // boundary_seq_start + index as u64,
                                    (alignment.xstart + index) as u64 + boundary_seq_start - 1,
                                    alignment.xstart + index,
                                    )
                                )
                            })
                            .filter(|(ena_id, pos, _, _, _, _, _)| {
                                !homology_loci
                                    .lock()
                                    .unwrap()
                                    .contains(&format!("{}_{}", ena_id, pos))
                            })
                            .collect::<Vec<_>>();

                        pos.iter()
                            .map(
                                |(ena_id, pos, s, ref_base, alt_base, genome_pos, xstart_index)| {
                                    homology_loci
                                        .lock()
                                        .unwrap()
                                        .insert(format!("{}_{}", ena_id, pos));
                                    (ena_id, pos, s, ref_base, alt_base, genome_pos, xstart_index)
                                },
                            )
                            .collect::<Vec<_>>();

                        let mut vcf_records: Vec<vcf::Record>;

                        // let vcf_records = pos
                        //     .iter()
                        //     .map(|(index, s, refbase, altbase, pos_assembly)| {});
                        if !pos.is_empty() {
                            let mut mut_items = pos.iter();
                            let mut last_item = mut_items.next().unwrap();
                            let mut res = vec![last_item.clone()];
                            let mut res_all = Vec::new();
                            if pos.len() == 1 {
                                res_all.push(res);
                            } else {
                                res = vec![last_item.clone()];
                                while let Some(current) = mut_items.next() {
                                    if current.1 == last_item.1 + 1 && current.2 == last_item.2 {
                                        res.push(current.clone());
                                    } else {
                                        if !res.is_empty() {
                                            res_all.push(res);
                                        }
                                        res = vec![current.clone()];
                                    }
                                    last_item = current;
                                }
                                if !res.is_empty() {
                                    res_all.push(res);
                                }
                            }

                            debug!("RES VARIANTS VECTOR = {:?}", &res_all);
                            vcf_records = res_all
                                .iter()
                                .map(|vars| {
                                    let initial_loci = vars[0].1;
                                    // let assembly_pos = (alignment.xstart + initial_loci.index) as i32;
                                    let assembly_pos = vars[0].5.clone() as i32;
                                    let refbase = vars[0].3.clone();
                                    let operation = &vars[0].2;
                                    let mut alternate_bases_missed: String = vars
                                        .iter()
                                        .flat_map(|(_, _, _, _, altbase, _, _)| altbase.chars())
                                        .collect();
                                    let mut alternate_bases = refbase.clone();
                                    alternate_bases_missed.push_str(&alternate_bases);
                                    alternate_bases = alternate_bases_missed;
                                    // alternate_bases.push_str(&alternate_bases_missed);

                                    let record = vcf::Record::builder()
                                        // .set_chromosome(chr_id.parse().unwrap())
                                        .set_chromosome(read_id.parse().unwrap())
                                        .set_position(Position::try_from(assembly_pos).unwrap())
                                        .set_reference_bases(refbase.parse().unwrap())
                                        .set_alternate_bases(alternate_bases.parse().unwrap())
                                        .build()
                                        .unwrap();

                                    record
                                    // Ok(())
                                })
                                // .map(|x| x)
                                .collect::<Vec<vcf::record::Record>>();

                            {
                                let mut vcf_records_to_write =
                                    vcf_records_to_write_child.lock().expect("Erro no vcf");
                                vcf_records_to_write.extend(vcf_records);
                            }

                            let boundary_cor_seqs = res_all
                                .iter()
                                .fold(Record::with_attrs(
                                    "initial",
                                    None,
                                    &boundary_seq,
                                ), |acc, vars| {
                                    let initial_loci = vars[0].1;
                                    // let assembly_pos = (alignment.xstart + initial_loci.index) as i32;
                                    let assembly_pos = vars[0].5.clone() as i32;
                                    let refbase = vars[0].3.clone();
                                    let operation = &vars[0].2;
                                    let mut alternate_bases_missed: String = vars
                                        .iter()
                                        .flat_map(|(_, _, _, _, altbase, _, xstart_index)| {
                                            altbase.chars()
                                        })
                                        .collect();
                                    // let mut alternate_bases = refbase;
                                    // alternate_bases.push_str(&alternate_bases_missed);
                                    let mut alternate_bases = refbase.clone();
                                    alternate_bases_missed.push_str(&alternate_bases);
                                    alternate_bases = alternate_bases_missed;
                                    let xstart_index = vars[0].6;
                                    let boundary_cor_seq = apply_variant(
                                        &acc.seq().to_vec(),
                                        xstart_index,
                                        alternate_bases.clone(),
                                        operation,
                                        is_revcomp,
                                    );
                                    // let read_id_cor = format!(
                                    //     "gbc{}_{}_var_{}|{}|{}",
                                    //     bf1.qseqid,
                                    //     bf2.qseqid,
                                    //     initial_loci,
                                    //     alternate_bases,
                                    //     vars[0].0,
                                    // );
                                    let mut split =vars[0].0.split("|");
                                    let read_id_cor = format!(
                                        "c_{}|{}|{}",
                                        initial_loci,
                                        alternate_bases,
                                        split.last().unwrap(),
                                    );
                                    let description_cor = None;
                                    // let sequence = b"ACGT";
                                    Record::with_attrs(
                                        &read_id_cor[..],
                                        description_cor,
                                        &boundary_cor_seq,
                                    )
                                });
                                // .collect::<Vec<Record>>();
                            
                            {
                                let mut output_fasta =
                                    output_fasta_child.lock().expect("Error on fasta");
                                let mut output_fasta_hash = output_fasta_hash_child
                                    .lock()
                                    .expect("Error on hash ena sequence id.");
                                
                                
                                if !output_fasta_hash.contains(&ena.id().to_string()) {
                                    let ena_id = ena.id().to_string();
                                    let ena_mod = Record::with_attrs(&ena_id[..], None, &ena.seq());
                                    output_fasta.push(ena_mod);
                                    
                                    output_fasta_hash.insert(ena_id);
                                }
                                // output_fasta.push(ena.clone());
                                // output_fasta.push(boundary_cor_seqs);
                                // output_fasta.push(record_bondary);
                                // if let Some(debug) = matches.value_of("debug") {
                                //     match hash_cds_fasta.get(&bf1.qseqid) {
                                //         Some(record) => output_fasta.push(record.clone()),
                                //         None => println!("Error, couldnt find sequence"),
                                //     }

                                //     match hash_cds_fasta.get(&bf2.qseqid) {
                                //         Some(record) => output_fasta.push(record.clone()),
                                //         None => println!("Error, couldnt find sequence"),
                                //     }
                                // };
                            }

                            let seqs_change = res_all
                                .iter()
                                .fold(BoundarySeq {
                                    contig: chr_id.parse().unwrap(),
                                    start_assembly: 0,
                                    end_assembly: 0,
                                    boundary_cor_seq: boundary_seq,
                                }, |acc, vars| {
                                    let initial_loci = vars[0].1;
                                    // let assembly_pos = (alignment.xstart + initial_loci.index) as i32;
                                    let assembly_pos = vars[0].5.clone() as i32;
                                    let refbase = vars[0].3.clone();
                                    let operation = &vars[0].2;
                                    let mut alternate_bases_missed: String = vars
                                        .iter()
                                        .flat_map(|(_, _, _, _, altbase, _, xstart_index)| {
                                            altbase.chars()
                                        })
                                        .collect();
                                    // let mut alternate_bases = refbase;
                                    // alternate_bases.push_str(&alternate_bases_missed);
                                    let mut alternate_bases = refbase.clone();
                                    alternate_bases_missed.push_str(&alternate_bases);
                                    alternate_bases = alternate_bases_missed;
                                    let xstart_index = vars[0].6;
                                    let boundary_seq = acc.boundary_cor_seq;
                                    let boundary_cor_seq = apply_variant(
                                        &boundary_seq,
                                        xstart_index,
                                        alternate_bases.clone(),
                                        operation,
                                        is_revcomp,
                                    );
                                    // let start_assembly = boundary_seq_start;
                                    // let end_assembly = boundary_seq_end;
                                    BoundarySeq {
                                        contig: chr_id.parse().unwrap(),
                                        start_assembly: boundary_seq_start.clone(),
                                        end_assembly: boundary_seq_end.clone(),
                                        boundary_cor_seq,
                                    }
                                });
                                // .collect::<Vec<BoundarySeq>>();
                            {
                                let mut seqs_to_change =
                                    seqs_to_change_child.lock().expect("Erro no fasta");
                                seqs_to_change.push(seqs_change);
                            }
                        }

                        let joined_genes =
                        blast_features.clone().into_iter().fold(HashSet::new(), |mut acc, (bf1, bf2)| {
                            acc.insert(bf1.qseqid);
                            acc.insert(bf2.qseqid);
                            acc});

                        warn!(
                            "Fixing INS/DEL between {:?} genes, in order to avoid disrupting frameshifts.", &joined_genes
                        )
                        // https://stackoverflow.com/questions/62253011/advance-next-to-peek-with-multipeek
                        // https://stackoverflow.com/questions/37684444/rust-iterators-and-looking-forward-peek-multipeek
                        // debug!("{:?}", &pos);
                    } else {
                        warn!(
                            "The pair {} and {} should remain as two separate CDS.",
                            &bf1.qseqid, &bf2.qseqid
                        );
                    }
                }

                // let _chunks_diff = dissimilar::diff(ena_seq_str, boundary_seq_str);

                // println!("CHUNKS = {:?}", &chunks_diff);
                // ... process the Uniprot entry ...
                // }
            } else {
                panic!("Error: protein subject is not from UniRef90 database.")
            }
        });

    dbg!("FASTA IS POISONED = ", output_fasta_ref.is_poisoned());
    dbg!("VCF IS POISONED = ", vcf_records_to_write_ref.is_poisoned());

    let output_fast_to_write = output_fasta_ref.lock().unwrap();

    info!(
        "Fixing {} possible frameshift errors between adjacent gene pairs.",
        &output_fast_to_write.len()
    );

    if matches.is_present("debug") {
        info!("Writing {} records to debug.fasta", &output_fast_to_write.len());
        let mut out_debug = PathBuf::from(&output_path_assembly_str);
        out_debug.set_file_name("debug.fasta");
        let debug_fasta = File::create(out_debug).unwrap();
        let mut writer = fasta::Writer::new(debug_fasta);

        for r in output_fast_to_write.iter() {
            writer
                .write_record(r)
                .expect("Error while writing FASTA output");
        }
        // .collect::<Vec<Result<vcf::record::Record, Box<dyn std::error::Error>>>>();
    };

    let mut vcf_records_to_write = vcf_records_to_write_ref.lock().unwrap();

    vcf_records_to_write.sort_by_key(|r: &vcf::Record| r.position());

    let mut header = vcf::Header::builder();
    for vcf_rec in vcf_records_to_write.iter() {
        let chr_id = vcf_rec.chromosome().to_string();
        if !contigs.contains(&chr_id) {
            header = header.add_contig(Contig::new(chr_id.clone()));
            // .build();

            contigs.insert(chr_id);
        }
    }

    // let header_final = header.build();
    // writer_vcf
    //     .write_header(&header_final)
    //     .expect("Coldnt write VCF file contig header.");

    // for vcf_rec in vcf_records_to_write.iter() {
    //     writer_vcf.write_record(&vcf_rec);
    // }

    let mut seqs_to_change = seqs_to_change_ref.lock().unwrap();
    seqs_to_change.par_sort_by(|a, b| b.start_assembly.cmp(&a.start_assembly));

    info!("Generating corrected Assembly.");
    let input_assembly = File::open(input_path_genome).unwrap();
    let assembly = fasta::Reader::new(&input_assembly);
    let mut new_fasta_assembly = Vec::new();
    let mut assembly_seq;
    for reco in assembly.records() {
        if let Ok(rec) = reco {
            assembly_seq = rec.seq().to_vec();
            for seq in seqs_to_change.iter() {
                if seq.contig == rec.id() {
                    assembly_seq = [
                        &assembly_seq[..seq.start_assembly as usize],
                        &seq.boundary_cor_seq,
                        &assembly_seq[seq.end_assembly as usize + 1..],
                    ]
                    .concat();
                }
            }
            let read_id_cor = format!("{}_cor", rec.id());
            let description_cor = None;
            // let sequence = b"ACGT";
            let new_contig_cor =
                Record::with_attrs(&read_id_cor[..], description_cor, &assembly_seq);

            new_fasta_assembly.push(new_contig_cor);
        }
    }

    
    let mut writer_cor = fasta::Writer::new(output_assembly);
    info!("Writing new assembly");
    for r in new_fasta_assembly.iter() {
        writer_cor
            .write_record(r)
            .expect("Error while writing FASTA output");
    }
    // dbg!(&seqs_to_change);

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
