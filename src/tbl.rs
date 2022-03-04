#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use bio::io::gff;
use multimap::MultiMap;
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{alpha1, digit1, line_ending, multispace1, newline, tab};
use nom::character::complete::{not_line_ending, space1};
use nom::character::streaming::char;
use nom::combinator::{map, map_res, not, opt, recognize};
use nom::multi::{fold_many0, many0, many0_count, many1, separated_list0};
use nom::number::complete::i32;
use nom::sequence::{preceded, separated_pair, terminated, tuple};
use nom::IResult;
use nom::Parser;
use nom::{error::ErrorKind, Err};
use nom_regex::str::re_find;
use nom_unicode::complete::alphanumeric1;
use noodles_gff::record::Strand;
use ordered_multimap::list_ordered_multimap::ListOrderedMultimap;
use std::fs::{self, File};
use std::io::{self, Write};

use anyhow::Context;
use std::fmt;
use std::path::Path;

/// A set of annotation table Records, each for chr/contig/sequence.
#[derive(Debug)]
pub struct Reader {
    records: Vec<Record>,
    n_records: usize,
}

/// An annotation table (.tbl) Record
#[derive(Debug)]
pub struct Record {
    pub seqid: String,
    pub table_id: String,
    pub features: Vec<Feature>,
}

/// An annotation table (.tbl) record
#[derive(Debug)]
pub struct Feature {
    pub start: u64,
    pub end: u64,
    pub feature_key: FeatureType,
    pub strand: Strand,
    pub qualifiers: ListOrderedMultimap<String, String>,
}

#[derive(Debug)]
pub enum FeatureType {
    CDS,
    Gene,
    TRNA,
    Undef,
}

impl fmt::Display for FeatureType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FeatureType::CDS => write!(f, "CDS"),
            FeatureType::Gene => write!(f, "gene"),
            FeatureType::TRNA => write!(f, "tRNA"),
            FeatureType::Undef => write!(f, "undefined"),
        }
    }
}

/// Utility function to nom parser that matches alphanumeric + underscore characters
fn alphanumeric_underscore(input: &str) -> IResult<&str, String> {
    let (input, parsed) = fold_many0(
        preceded(opt(char('_')), alphanumeric1),
        Vec::new,
        |mut acc: Vec<_>, nxt: &str| {
            acc.push(nxt);
            acc
        },
    )(input)?;
    Ok((input, parsed.join("_")))
}

/// Utility function to nom parser that matches a string of characters up to any whitespace
fn everything_but_whitespace(input: &str) -> IResult<&str, String> {
    let re = regex::Regex::new(r"[^\s]+").unwrap();
    let parser = re_find::<(&str, ErrorKind)>(re);
    let (input, parsed) = parser(input).unwrap();
    Ok((input, parsed.to_string()))
}

// fn alphanumeric_doispontos(input: &str) -> IResult<&str, String> {
//     let (input, parsed) = fold_many0(
//         preceded(opt(char(':')), alphanumeric1),
//         Vec::new,
//         |mut acc: Vec<_>, nxt: &str| {
//             acc.push(nxt);
//             acc
//         },
//     )(input)?;
//     Ok((input, parsed.join(":")))
// }

fn only_seqid(input: &str) -> IResult<&str, (String, String)> {
    let (input, parsed) = fold_many0(
        preceded(opt(char('_')), alphanumeric1),
        Vec::new,
        |mut acc: Vec<_>, nxt: &str| {
            acc.push(nxt);
            acc
        },
    )(input)?;
    Ok((input, (parsed.join("_"), "".to_string())))
}

/// Parse into a TBL Record struct.
pub fn parse_input(input: &str) -> IResult<&str, Record> {
    let (input, _) = tag(">Feature ")(input)?;
    // let initial = tag(">Feature ");
    let both = separated_pair(everything_but_whitespace, space1, everything_but_whitespace);
    let only_one = alphanumeric_underscore;
    let (input, (seqid, table_name)) = terminated(alt((both, only_seqid)), newline)(input)?;

    let feat = tuple((
        terminated(digit1, tab),
        terminated(digit1, tab),
        terminated(alpha1, newline),
    ));
    let qualifier = tuple((tab, tab, tab, alphanumeric_underscore, tab, not_line_ending));
    let qualifiers = map(separated_list0(newline, qualifier), |x| {
        let map_features: ListOrderedMultimap<String, String> = x
            .iter()
            .map(|(_, _, _, qualifier_key, _, qualifier_value)| {
                (qualifier_key.to_string(), qualifier_value.to_string())
            })
            .collect();
        map_features
    });

    let snippet = map(
        tuple((feat, qualifiers)),
        |((start, end, feat_type), vec_qualifiers)| {
            let pseudo_start: u64 = start.parse().unwrap();
            let pseudo_end: u64 = end.parse().unwrap();
            let feat_type_enum = match feat_type {
                "CDS" => FeatureType::CDS,
                "gene" => FeatureType::Gene,
                "tRNA" => FeatureType::TRNA,
                _ => FeatureType::Undef,
            };

            // Features that are on the complementary strand are indicated by reversing the invertal locations
            if pseudo_start > pseudo_end {
                Feature {
                    start: end.parse().unwrap(),
                    end: start.parse().unwrap(),
                    strand: Strand::Reverse,
                    feature_key: feat_type_enum,
                    qualifiers: vec_qualifiers,
                }
            } else {
                Feature {
                    start: start.parse().unwrap(),
                    end: end.parse().unwrap(),
                    strand: Strand::Forward,
                    feature_key: feat_type_enum,
                    qualifiers: vec_qualifiers,
                }
            }
        },
    );

    let (input, parsed_feats) = many1(terminated(snippet, opt(line_ending)))(input)?;

    // dbg!(parsed_feats);

    Ok((
        input,
        Record {
            seqid: seqid,
            table_id: table_name,
            features: parsed_feats,
        },
    ))
}

/// Parse a set of Records into a TBL Reader struct.
/// A TBL may contain a set of annotations (features) for each sequence/contig/chromossome.
pub fn parse_file(input: &str) -> IResult<&str, Reader> {
    let (input, parsed) = many1(parse_input)(input)?;
    Ok((
        input,
        Reader {
            n_records: parsed.len(),
            records: parsed,
        },
    ))
}

impl Record {
    pub fn new<R: std::io::Read>(reader: &mut R) -> Self {
        let mut tbl_string = String::new();
        reader
            .read_to_string(&mut tbl_string)
            .expect("error while parsing record");
        let (_, parsed) = parse_input(&tbl_string[..]).unwrap();
        parsed
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Feature> {
        self.features.iter_mut()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Feature> {
        self.features.iter()
    }
}

impl Reader {
    /// Create a new TBL reader given an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # extern crate rs_genomesub;
    /// # use crate::rs_genomesub::tbl::Reader;
    /// # fn main() {
    /// const TBL_FILE: &'static [u8] = b">Feature Ps_genome
    /// 219	617	CDS
    /// db_xref	COG:COG0012
    /// gene	ychF_1
    /// inference	ab initio prediction:Prodigal:002006";
    /// let reader = Reader::new(&mut TBL_FILE);
    /// # }
    /// ```
    pub fn new<R: std::io::Read>(reader: &mut R) -> Self {
        let mut tbl_string = String::new();
        reader
            .read_to_string(&mut tbl_string)
            .expect("Error while creating new record");
        let (_, parsed) = parse_file(&tbl_string[..]).unwrap();
        parsed
    }

    /// Read TBL from given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(|mut f| Reader::new(&mut f))
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Record> {
        self.records.iter_mut()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Record> {
        self.records.iter()
    }
}

// pub trait ToFile {
//     /// a mutable reference to any value that implements the Write trait
//     fn send_to_file<W: Write>(&self, out: &mut W) -> io::Result<()>;

//     fn convert_to_gff<W: Write>(&self, writer: W) -> io::Result<()>;
// }

// impl ToFile for Record {
//     fn send_to_file<W: Write>(&self, out: &mut W) -> io::Result<()> {
//         if self.table_id.is_empty() {
//             writeln!(out, ">Feature {}", self.seqid)?;
//         } else {
//             writeln!(out, ">Feature {} {}", self.seqid, self.table_id)?;
//         }

//         for feat in &self.features {
//             match feat.strand {
//                 Strand::Forward => {
//                     writeln!(out, "{}\t{}\t{}", feat.start, feat.end, feat.feature_key)?;
//                 }

//                 Strand::Reverse => {
//                     writeln!(out, "{}\t{}\t{}", feat.end, feat.start, feat.feature_key)?;
//                 }

//                 Strand::Unknown => {
//                     panic!("Could not identify strand of feature {:?}", &feat);
//                 }
//             }

//             // feat.qualifiers.keys().collect::<Vec<String>>();

//             // for key in feat.qualifiers.keys() {
//             //     if feat.qualifiers.is_vec(key) {
//             //         for q in feat.qualifiers.get_vec(key).unwrap() {
//             //             writeln!(out, "\t\t\t{}\t{}", key, q)?;
//             //         }
//             //     } else {
//             //         writeln!(out, "\t\t\t{}\t{}", key, feat.qualifiers[key])?;
//             //     }
//             // }
//             for (qual_key, qual_val) in feat.qualifiers.iter() {
//                 writeln!(out, "\t\t\t{}\t{}", qual_key, qual_val)?;
//             }
//         }

//         Ok(())
//     }

//     fn convert_to_gff<W: Write>(&self, writer: W) -> io::Result<()> {
//         // let mut records: Vec<Record> = Vec::new();

//         let mut writergff = gff::Writer::new(writer, gff::GffType::GFF3);

//         for feat in &self.features {
//             let mut gffrecord = gff::Record::new();

//             let mut seqname = gffrecord.seqname_mut();
//             *seqname = self.seqid.clone();

//             let mut source = gffrecord.source_mut();
//             *source = ".".to_string();

//             let mut feature_type = gffrecord.feature_type_mut();
//             *feature_type = feat.feature_key.clone();

//             let mut start = gffrecord.start_mut();
//             *start = feat.start;

//             let mut end = gffrecord.end_mut();
//             *end = feat.end;

//             let mut score = gffrecord.score_mut();
//             *score = ".".to_string();

//             let mut strand = gffrecord.strand_mut();
//             *strand = (feat.strand.strand_symbol()).to_string();

//             let mut frame = gffrecord.frame_mut();
//             *frame = ".".to_string();

//             let mut attributes = MultiMap::new();
//             for (k, v) in &feat.qualifiers {
//                 attributes.insert(k.clone(), v.clone());
//             }

//             let mut att = gffrecord.attributes_mut();
//             *att = attributes;

//             // records.push(gffrecord);
//             writergff
//                 .write(&gffrecord)
//                 .ok()
//                 .expect("Error writing record.");
//         }

//         Ok(())
//     }
// }

#[cfg(test)]
mod tests {

    use super::*;

    const TBL_COMPACT: &[u8] = b">Feature Ps_genome
219	617	CDS
			db_xref	COG:COG0012
			gene	ychF_1
			inference	ab initio prediction:Prodigal:002006
			inference	similar to AA sequence:UniProtKB:P0ABU2
			locus_tag	LZT29_00001
			product	Ribosome-binding ATPase YchF
827	1264	CDS
			db_xref	COG:COG0582
			gene	intS_1
			inference	ab initio prediction:Prodigal:002006
			inference	similar to AA sequence:UniProtKB:P37326
			locus_tag	LZT29_00002
			product	Prophage integrase IntS
";

    #[test]
    fn test_file_open() {
        let path = "examples/input_multi.tbl";

        let res = Reader::from_file(path);

        assert_eq!(res.unwrap().n_records, 2);
    }

    #[test]
    fn test_reader() {
        // let path = "examples/input_text.txt";
        // let mut input = File::open(path).unwrap();

        // let data = b">Feature Ps_genome tabelax\n";
        let mut reader = &TBL_COMPACT[..];

        let res = Record::new(&mut reader);

        // let mut handler = File::create("tmp.txt").unwrap();
        // res.send_to_file(&mut handler);

        assert_eq!(res.features.len(), 2);
        assert_eq!(
            res.features[0].qualifiers.get(&"locus_tag".to_string()),
            Some(&"LZT29_00001".to_string())
        );
        assert_eq!(
            res.features[1].qualifiers.get(&"locus_tag".to_string()),
            Some(&"LZT29_00002".to_string())
        );

        // for feat in &res.features {
        //     println!(
        //         "GO {:?}",
        //         feat.qualifiers.get(&"locus_tag".to_string()).unwrap()
        //     )
        // }

        // println!("WHYY GEORGIA!");
    }
}
