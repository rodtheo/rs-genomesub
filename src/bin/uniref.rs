extern crate libflate;
extern crate uniprot;
extern crate ureq;

use bio::io::fasta::{self, IndexedReader, Record};
// use nom::error::dbg_dmp;
// use rs_genomesub::tbl::Feature;
// use uniprot::uniprot::DbReference;
use ureq::Response;

extern crate rs_genomesub;
use crate::rs_genomesub::blast_tabular::*;
use crate::rs_genomesub::tbl;

use itertools::Itertools;
use std::collections::HashMap;

use std::fs::File;
// use std::io::Write;
use std::path::Path;
// use std::path::PathBuf;

use protein_translate::translate;

// Record name
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

pub fn get_record_from_insilico_translation(
    response: Response,
    prot_id: &str,
) -> Result<Record, std::io::Error> {
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

    println!("STRING SEQ: {}", &protseq);

    // println!("{:?}", &entry);

    // let protseq = b"ATC";
    // let read_id = "ehissoai";
    let read_id = format!("InSilico_{}", prot_id);
    let protseq_u8 = protseq.as_bytes();
    let dna = translate(protseq_u8);

    dbg!("UAI, ", &dna);

    let description = None;
    let record = Record::with_attrs(&read_id, description, dna.as_bytes());

    dbg!("{:?}", &record);

    Ok(record)
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

    let mut input = File::open("genes_out.fa").unwrap();
    let res = fasta::Reader::new(&mut input);

    let mut hash_cds_fasta = HashMap::new();

    let mut record;
    let mut record_id;
    for s in res.records() {
        record = s.unwrap();
        record_id = record.id().to_string();
        hash_cds_fasta.insert(record_id, record);
    }

    let mut input_tbl = File::open("pantoea_genome_add_gene_feats.tbl").unwrap();
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

    let mut output_fasta: Vec<fasta::Record> = Vec::new();

    let fblastout = include_str!("../../examples/output_diamond_sample.txt");

    let (_inp, res) = parse_input(fblastout).unwrap();

    let path_fsa = Path::new("pantoea_genome_add_gene_feats.fsa");

    let mut faidx = IndexedReader::from_file(&path_fsa).unwrap();

    // let count_possible_frameshifts = res
    //     .into_iter()
    //     .tuple_windows()
    //     .filter(|(s1, s2)| s1.sseqid == s2.sseqid)
    //     .count();

    let count_possible_frameshifts = res
        .into_iter()
        .tuple_windows()
        .filter(|(s1, s2)| (s1.sseqid == s2.sseqid) & (s1.qframe == s2.qframe))
        .collect::<Vec<(BlastFeature, BlastFeature)>>();

    eprintln!(
        "Number of pairs with different frameshift matching the same gene: {}",
        &count_possible_frameshifts.len()
    );

    for (bf1, bf2) in count_possible_frameshifts.iter() {
        if let Some(query) = bf1.sseqid.strip_prefix("UniRef90_") {
            let query_url = format!(
                "https://www.uniprot.org/uniprot/{}.xml?include=yes&compress=yes",
                // "https://www.uniprot.org/uniref/UniRef90_{}.xml?include=yes&compress=yes",
                query
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

                    eprintln!(
                        "Error when mapping to reciprocal ProteinKB entry for UniRef90_{}",
                        &query
                    );
                    get_record_from_insilico_translation(req, query)
                }
            };
            // let reader = libflate::gzip::Decoder::new(req.call().unwrap().into_reader()).unwrap();

            // for r in uniprot::uniprot::parse(std::io::BufReader::new(reader)) {

            match hash_cds_fasta.get(&bf1.qseqid) {
                Some(record) => output_fasta.push(record.clone()),
                None => println!("Error, couldt find sequence"),
            }

            match hash_cds_fasta.get(&bf2.qseqid) {
                Some(record) => output_fasta.push(record.clone()),
                None => println!("Error, couldt find sequence"),
            }

            // let ena_record = reader.records().next().unwrap();

            let ena = match reader {
                Ok(rec) => {
                    output_fasta.push(rec.clone());
                    rec
                }
                Err(_e) => panic!("Could find Record"),
            };

            let ena_seq_str = std::str::from_utf8(ena.seq()).unwrap();

            let (tbl_bf1, chr_id) = hash_tbl_feat.get(&bf1.qseqid).unwrap();
            let (tbl_bf2, chr_id2) = hash_tbl_feat.get(&bf2.qseqid).unwrap();

            assert_eq!(chr_id, chr_id2);
            let boundary_seq_description;
            let boundary_seq = if tbl_bf1.start < tbl_bf2.start {
                boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf1.start, tbl_bf2.end);
                get_genomic_sequence(chr_id, tbl_bf1.start, tbl_bf2.end, &mut faidx)
            } else {
                boundary_seq_description = format!("{}:{}-{}", chr_id, tbl_bf2.start, tbl_bf1.end);
                get_genomic_sequence(chr_id, tbl_bf2.start, tbl_bf1.end, &mut faidx)
            };

            let read_id = format!("my_genomic_boundary_{}_{}", bf1.qseqid, bf2.qseqid);
            let description = Some(&boundary_seq_description[..]);
            // let sequence = b"ACGT";
            let record = Record::with_attrs(&read_id[..], description, &boundary_seq);
            output_fasta.push(record);

            let boundary_seq_str = std::str::from_utf8(&boundary_seq).unwrap();

            let _chunks_diff = dissimilar::diff(ena_seq_str, boundary_seq_str);

            // println!("CHUNKS = {:?}", &chunks_diff);
            // ... process the Uniprot entry ...
            // }
        } else {
            panic!("Error: protein subject is not from UniRef90 database.")
        }
    }

    let mut writer = fasta::Writer::new(std::io::stdout());

    for r in output_fasta.iter() {
        writer
            .write_record(r)
            .expect("Error while writing FASTA output");
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
