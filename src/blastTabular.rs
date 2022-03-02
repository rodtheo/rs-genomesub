use nom::combinator::{map, opt};
use nom::error::ErrorKind;
use nom::IResult;
use nom_regex::str::re_find;

use nom::branch::alt;
use nom::bytes::complete::take;
use nom::character::complete::{line_ending, multispace0, multispace1, newline, tab};
use nom::error::Error;
use nom::multi::{many1, separated_list1};
use nom::number::complete::{u16, u64};
use nom::sequence::{terminated, tuple};

use itertools::Itertools;

#[derive(Debug, Default, Clone)]
pub struct BlastFeature {
    // qseqid Query Seq - id
    pub qseqid: String,
    // stitle Subject Title
    pub stitle: String,
    // sseqid Subject Seq - id
    pub sseqid: String,
    // qstart Start of alignment in query*
    pub qstart: u64,
    // qend End of alignment in query*
    pub qend: u64,
    // sstart Start of alignment in subject*
    pub sstart: u64,
    // send End of alignment in subject*
    pub send: u64,
    // qframe Query frame
    pub qframe: i16,
    // btop Blast traceback operations(BTOP)*
    pub btop: String,
}

/// Utility function to nom parser that matches a string of characters up to a tab or new line
fn everything_but_tab_or_nl(input: &str) -> IResult<&str, String> {
    let re = regex::Regex::new(r"[^\t\n]+").unwrap();
    let parser = re_find::<(&str, ErrorKind)>(re);
    let (input, parsed) = parser(input).unwrap_or((input, ""));
    Ok((input, parsed.to_string()))
}

// fn parse_u64(input: &str) -> IResult<&str, String> {
// 	let parser =
// }

pub fn parse_input(input: &str) -> IResult<&str, Vec<BlastFeature>> {
    let line = tuple((
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        tab,
        everything_but_tab_or_nl,
        line_ending,
    ));

    let map_qualifiers = map(line, |x| BlastFeature {
        qseqid: x.0,
        stitle: x.2,
        sseqid: x.4,
        qstart: x.6.parse::<u64>().unwrap(),
        qend: x.8.parse::<u64>().unwrap(),
        sstart: x.10.parse::<u64>().unwrap(),
        send: x.12.parse::<u64>().unwrap(),
        qframe: x.14.parse::<i16>().unwrap(),
        btop: x.16,
    });

    // let term = terminated(map_qualifiers, opt(line_ending));

    // dbg!(map_qual_line);

    // let feats = separated_list1(line_ending, map_qualifiers);

    let feats = many1(map_qualifiers);

    let (input, qualifiers) = terminated(feats, opt(multispace1))(input)?;

    // dbg!(&qualifiers);

    // let (input, qualifiers) = map(separated_list1(newline, line), |x| {
    //     let map_features: Vec<BlastFeature> = x
    //         .iter()
    //         .map(
    //             |(qseqid, _, stitle, _, _, _, _, _, _, _, _, _, _, _, _, _, _)| {
    //                 BlastFeature::default()
    //             },
    //         )
    //         .collect::<Vec<BlastFeature>>();
    //     map_features
    // })(input)?;

    // let (input, qualifiers) = many1(line)(input)?;

    // dbg!(input);

    Ok((input, qualifiers))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_open() {
        let fblastout = include_str!("../out_diamond.txt");

        let (inp, res) = parse_input(fblastout).unwrap();

        dbg!(res);

        // let res = Reader::from_file(path);
        assert_eq!(1, 1);
        // assert_eq!(res.unwrap().n_records, 2);
    }

    #[test]
    fn test_iter() {
        let fblastout = include_str!("../examples/output_diamond_sample.txt");

        let (inp, res) = parse_input(fblastout).unwrap();

        // let count_possible_frameshifts = res
        //     .into_iter()
        //     .tuple_windows()
        //     .filter(|(s1, s2)| s1.sseqid == s2.sseqid)
        //     .count();

        let count_possible_frameshifts = res
            .into_iter()
            .tuple_windows()
            .filter(|(s1, s2)| (s1.sseqid == s2.sseqid) & (s2.qframe == 2))
            .collect::<Vec<(BlastFeature, BlastFeature)>>();

        dbg!(count_possible_frameshifts);

        assert_eq!(1, 1);
    }
}
