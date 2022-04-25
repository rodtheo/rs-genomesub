use bio::alphabets::dna;
use noodles::gff::record::Strand;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use noodles::gff::{
    self as gff,
    record::{attributes::Entry, Attributes},
};

use bio::io::fasta::{self, IndexedReader};
use std::path::Path;

use clap::{arg, Arg, Command};

extern crate rs_genomesub;
use crate::rs_genomesub::tbl::*;

/// Tool for processing NCBI .tbl file format.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let input = include_str!("../input_text.txt").trim();

    let matches = Command::new("tbltk")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .about(
            "A toolkit to manipulate NCBI's Feature Table Format (suffix .tbl).\
            \nSpecs on https://www.ncbi.nlm.nih.gov/genbank/feature_table/",
        )
        .version("v0.1.0")
        // .author("Rodrigo Theodoro Rocha <theodorobiotec@gmail.com>")
        // .subcommand(subcmd)
        .subcommand(
            Command::new("togff")
                .arg_required_else_help(true)
                .about("Converts a TBL to a GFF3")
                .arg(arg!(-g --genome <GENOME>).required(true))
                .arg(arg!([input] "TBL file to operate on").required(true)),
        )
        .subcommand(
            Command::new("tofasta")
                .arg_required_else_help(true)
                .about("Converts a tbl annotations into nucleotide gene sequences in fasta")
                .arg(arg!(-g --genome <GENOME>).required(true))
                .arg(arg!([input] "TBL file to operate on").required(true)),
        )
        // .arg(
        //     // Arg::with_name("input")
        //     //     .value_name("INPUT")
        //     //     .help("Input tbl")
        //     arg!([input] "TBL file to operate on").required(true), // .min_values(1),
        // )
        // .arg(
        //     Arg::new("genome")
        //         .value_name("GENOME")
        //         .help("Genome in FASTA. Required if --to-fast is choose.")
        //         .takes_value(true),
        // )
        // // .arg(
        // //     Arg::new("togff")
        // //         .value_name("TOGFF")
        // //         .help("Converts a TBL to a GFF3")
        // //         .takes_value(false)
        // //         .short('g')
        // //         .long("to-gff")
        // //         .conflicts_with("tofasta"),
        // // )
        // .arg(
        //     Arg::new("tofasta")
        //         .value_name("TOFASTA")
        //         .help("Converts a TBL to a nucleotide gene sequences in FASTA format")
        //         .takes_value(false)
        //         .short('f')
        //         .long("to-fasta")
        //         .conflicts_with("togff")
        //         .requires("genome"),
        // )
        .get_matches();

    // let input_path_str = matches.value_of("input").unwrap();
    // let togff = matches.is_present("togff");
    // let tofasta = matches.is_present("tofasta");

    // let input_path = PathBuf::from(&input_path_str);
    // let mut output_path = PathBuf::from(&input_path_str);

    match matches.subcommand() {
        Some(("tofasta", sub_matches)) => {
            let input_path_str = sub_matches.value_of("input").unwrap();
            let input_path = PathBuf::from(&input_path_str);

            let mut input = File::open(input_path).unwrap();
            let res = Reader::new(&mut input);

            let mut output_path = PathBuf::from(&input_path_str);
            output_path.set_extension("fa");

            let genome_path_str = sub_matches.value_of("genome").unwrap();

            let path_fsa = Path::new(genome_path_str);

            // let mut faidx = IndexedReader::from_file(&path_fsa).expect(&format!("Please, provide an index for {}. It can be accomplished by running `samtools faidx {}`", &genome_path_str, &genome_path_str));
            let faidx = IndexedReader::from_file(&path_fsa);
            let mut faidx = match faidx {
            Ok(file) => file,
            Err(_) => return Err(format!("Couldn't open the sequence file index! Please, provide an index for {}. It can be created by running `samtools faidx {}`", &genome_path_str, &genome_path_str).into()),
            };

            let mut seq = Vec::new();

            let mut writer = fasta::Writer::new(std::io::stdout());

            // faidx.read(&mut seq).expect("Coldn't read the interval");
            // println!("{:?}", &seq);

            // faidx
            //     .fetch("Ps_genome", 100, 110)
            //     .expect("Couldn't fetch interval");
            // faidx.read(&mut seq).expect("Coldn't read the interval");
            // println!("{:?}", &seq);
            for record in res.iter() {
                let count = record
                    .iter()
                    .filter(|f| matches!(f.feature_key, FeatureType::CDS))
                    .count();
                eprintln!("Number of CDS = {}", count);
            }

            for record in res.iter() {
                for feat in record
                    .iter()
                    .filter(|f| matches!(f.feature_key, FeatureType::CDS))
                {
                    faidx
                        .fetch(&record.seqid[..], feat.start - 1, feat.end)
                        .expect("Couldn't fetch interval");
                    faidx.read(&mut seq).expect("Couldn't read the interval");

                    if matches!(feat.strand, Strand::Reverse) {
                        seq = dna::revcomp(seq);
                    }

                    writer
                        .write(
                            feat.qualifiers.get("locus_tag").unwrap(),
                            None,
                            seq.as_slice(),
                        )
                        .expect("Error writing fasta record.");
                    // println!("{:?}", &feat);
                    // if feat.feature_key
                    // faidx
                    //     .fetch("Ps_genome", 0, 100)
                    //     .expect("Couldn't fetch interval");
                    // // println!("{:?}", &feat)
                }
            }
        }

        Some(("togff", sub_matches)) => {
            let input_path_str = sub_matches.value_of("input").unwrap();
            let input_path = PathBuf::from(&input_path_str);

            let mut input = File::open(input_path).unwrap();
            let res = Reader::new(&mut input);

            let mut output_path = PathBuf::from(&input_path_str);
            output_path.set_extension("gff");

            let mut writer = gff::Writer::new(Vec::new());
            let version = gff::Directive::GffVersion(Default::default());
            writer.write_directive(&version)?;

            for record in res.iter() {
                for feat in record.iter() {
                    let mut vec_attributes: Vec<Entry> = Vec::new();
                    for (k, v) in &feat.qualifiers {
                        vec_attributes.push(Entry::new(k.clone(), v.clone()));
                    }

                    let gff_record = gff::Record::builder()
                        .set_reference_sequence_name(record.seqid.clone())
                        .set_start(feat.start as i32)
                        .set_end(feat.end as i32)
                        .set_type(feat.feature_key.to_string())
                        .set_strand(feat.strand)
                        .set_attributes(Attributes::from(vec_attributes))
                        .build();

                    writer
                        .write_record(&gff_record)
                        .expect("Error writing gff record");
                    // println! {"{:?}\n\n", &feat}
                }
            }

            // let record = gff::Record::builder()
            //                         .set_reference_sequence_name(String::from("sq0"))
            //                         .build();

            // let record = gff::Record::default();
            // writer.write_record(&record)?;

            let mut file = File::create(output_path)?;
            // let writer_file = BufWriter::new(file);
            file.write_all(writer.get_ref())
                .expect("Unable to write data.");
        }

        _ => unreachable!("Exhausted list of subcommands and subcommand_required prevents `None`"),
    }

    Ok(())
}
