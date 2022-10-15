#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]

use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;
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
use noodles::gff::record::Strand;
use ordered_multimap::list_ordered_multimap::ListOrderedMultimap;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Write};

use anyhow::Context;
use std::fmt;
use std::path::Path;

use flate2::read::GzDecoder;
use std::io::prelude::*;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;

use super::*;

/// A set of annotation table Records, each for chr/contig/sequence.
#[derive(Debug)]
pub struct Chains {
    pub records: Vec<Chain>,
    pub n_records: usize,
}

impl Chains {
    pub fn to_lift(&self) -> AnnotMap<String, liftover::LiftBlock> {
        // let mut chains: HashMap<String, liftover::LiftBlock> = HashMap::new();
        let mut chains = AnnotMap::new();

        for chain in self.records.iter() {
            let coords_ins = chain
                .alignments
                .iter()
                .flat_map(|aln| {
                    if let (Some(dt), Some(dq)) = (aln.dt, aln.dq) {
                        // Some(aln.size + dq)
                        std::iter::repeat(Match)
                            .take(aln.size)
                            .chain(std::iter::repeat(Del).take(dt))
                            .chain(std::iter::repeat(Ins).take(dq))
                    } else {
                        std::iter::repeat(Ins)
                            .take(aln.size)
                            .chain(std::iter::repeat(Del).take(0))
                            .chain(std::iter::repeat(Ins).take(0))
                    }
                })
                .collect::<Vec<AlignmentOperation>>();

            // dbg!(&coords_ins.len());

            let strand = if chain.tStrand == '+' {
                ReqStrand::Forward
            } else {
                ReqStrand::Reverse
            };

            let contig_seq_entry = Contig::new(
                chain.tName.clone(),
                chain.tStart as isize,
                (chain.tEnd - chain.tStart),
                strand,
            );

            let block = liftover::LiftBlock::new(&coords_ins);

            chains.insert_at(block, &contig_seq_entry);
            // chains.insert(chain.tName.clone(), block);
        }

        chains
    }
}

/// An annotation Chain
#[derive(Debug)]
pub struct Chain {
    // The alignment start and end positions are represented as zero-based half-open intervals.
    // For example, the first 100 bases of a sequence would be represented with start position = 0
    // and end position = 100, and the next 100 bases would be represented as start position = 100
    //  and end position = 200. When the strand value is "-", position coordinates are listed in
    // terms of the reverse-complemented sequence.

    // chain score
    pub score: usize,
    // chromosome (reference/target sequence)
    pub tName: String,
    // chromosome size (reference/target sequence)
    pub tSize: usize,
    // strand (reference/target sequence)
    pub tStrand: char,
    // alignment start position (reference/target sequence)
    pub tStart: usize,
    // alignment end position (reference/target sequence)
    pub tEnd: usize,
    // chromosome (query sequence)
    pub qName: String,
    // chromosome size (query sequence)
    pub qSize: usize,
    // strand (query sequence)
    pub qStrand: char,
    // alignment start position (query sequence)
    pub qStart: usize,
    // alignment end position (query sequence)
    pub qEnd: usize,
    // chain ID
    pub id: usize,
    // Alignment Data Lines
    pub alignments: Vec<ChainAlignment>,
}

///  Alignment data lines contain three required attribute values: size, dt and dq
//  NOTE: The last line of the alignment section contains only one number:
// the ungapped alignment size of the last block.
#[derive(Debug)]
pub struct ChainAlignment {
    // the size of the ungapped alignment
    pub size: usize,
    // the difference between the end of this
    // block and the beginning of the next block
    //  (reference/target sequence)
    pub dt: Option<usize>,
    // the difference between the end of this block
    //  and the beginning of the next block (query sequence)
    pub dq: Option<usize>,
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

/// Parse into a Chain Record struct.
pub fn parse_input(input: &str) -> IResult<&str, Chain> {
    let (input, _) = tag("chain ")(input)?;
    let mut header_line = tuple((
        terminated(digit1, space1),
        terminated(everything_but_whitespace, space1),
        terminated(digit1, space1),
        terminated(alt((tag("+"), tag("-"))), space1),
        terminated(digit1, space1),
        terminated(digit1, space1),
        terminated(everything_but_whitespace, space1),
        terminated(digit1, space1),
        terminated(alt((tag("+"), tag("-"))), space1),
        terminated(digit1, space1),
        terminated(digit1, space1),
        terminated(digit1, newline),
    ));

    let (input, hl) = header_line(input)?;

    let alignment_line = tuple((
        terminated(digit1, opt(space1)),
        opt(terminated(digit1, space1)),
        opt(digit1),
    ));

    let mut alignment_lines = map(
        separated_list0(newline, alignment_line),
        |x: Vec<(&str, Option<&str>, Option<&str>)>| {
            x.iter()
                .map(|(size, dt, dq)| ChainAlignment {
                    size: size.parse().unwrap(),
                    dt: match dt {
                        Some(val) => Some(val.parse::<usize>().unwrap()),
                        None => None,
                    },
                    dq: match dq {
                        Some(val) => Some(val.parse::<usize>().unwrap()),
                        None => None,
                    },
                })
                .collect::<Vec<ChainAlignment>>()
        },
    );

    let (input, al) = alignment_lines(input)?;

    let (input, _) = newline(input)?;

    let aln_chain = Chain {
        score: hl.0.parse().unwrap(),
        tName: hl.1.parse().unwrap(),
        tSize: hl.2.parse().unwrap(),
        tStrand: hl.3.parse().unwrap(),
        tStart: hl.4.parse().unwrap(),
        tEnd: hl.5.parse().unwrap(),
        qName: hl.6.parse().unwrap(),
        qSize: hl.7.parse().unwrap(),
        qStrand: hl.8.parse().unwrap(),
        qStart: hl.9.parse().unwrap(),
        qEnd: hl.10.parse().unwrap(),
        id: hl.11.parse().unwrap(),
        alignments: al,
    };

    Ok((input, aln_chain))
}

pub fn parse_file(input: &str) -> IResult<&str, Chains> {
    let (input, parsed) = many1(terminated(parse_input, opt(newline)))(input)?;
    let n_entries = parsed.len();
    Ok((
        input,
        Chains {
            records: parsed,
            n_records: n_entries,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_open() {
        let chainout = include_str!("../out.chain");

        let (inp, res) = parse_input(chainout).unwrap();

        // dbg!(inp);
        // dbg!(res);

        // let res = Reader::from_file(path);
        assert_eq!(1, 1);
        // assert_eq!(res.unwrap().n_records, 2);
    }

    #[test]
    fn test_file_multiple_chains_open() {
        let p = Path::new("examples/hg19ToHg38.over.chain.gz");
        // dbg!(p);
        let mut reader = File::open(p).map(GzDecoder::new).unwrap();

        let mut data = String::new();
        reader.read_to_string(&mut data).unwrap();

        let (inp, res) = parse_file(&data[..]).unwrap();

        assert_eq!(res.n_records, 1278);
    }
}
