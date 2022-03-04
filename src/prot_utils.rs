use multimap::MultiMap;

use nom::combinator::{map, opt};
use nom::error::ErrorKind;
use nom::IResult;
use nom_regex::str::re_find;

use nom::bytes::complete::{tag, take_till};
use nom::character::complete::{digit1, line_ending, newline, not_line_ending, space0, space1};
use nom::multi::{many1, many_till, separated_list1};
use nom::sequence::{terminated, tuple};

#[derive(Debug, Default, Clone)]
pub struct GeneticCodeTable<'a> {
    name: String,
    shortname: Option<String>,
    id: u64,
    start_codons: Vec<&'a str>,
    end_codons: Vec<&'a str>,
    table: MultiMap<u8, Vec<u8>>,
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

fn till_aspas(s: &str) -> IResult<&str, &str> {
    take_till(|c| c == '"')(s)
}

pub fn parse_input(input: &str) -> IResult<&str, Vec<GeneticCodeTable>> {
    let comments = tuple((tag("--"), not_line_ending, line_ending));

    let (input, (_, _line_rest)) = many_till(comments, tag("\nGenetic-code-table ::= {\n"))(input)?;

    let begin_bracket = tuple((space1, tag("{"), space0, newline));

    let table_line_name = tuple((
        space1,
        tag("name"),
        space1,
        tag("\""),
        till_aspas,
        tag("\""),
        space1,
        tag(","),
        line_ending,
    ));

    let table_line_names = many1(table_line_name);
    let table_line_id = tuple((space1, tag("id"), space1, digit1, space1, tag(",")));
    let table_line_prot = tuple((
        space1,
        tag("ncbieaa"),
        space1,
        tag("\""),
        till_aspas,
        tag("\""),
        space0,
        tag(","),
        newline,
    ));
    let table_line_s = tuple((space1, tag("sncbieaa"), space1, not_line_ending));
    let table_line_base1 = tuple((space1, tag("-- Base1"), space1, not_line_ending));
    let table_line_base2 = tuple((space1, tag("-- Base2"), space1, not_line_ending));
    let table_line_base3 = tuple((space1, tag("-- Base3"), space1, not_line_ending));
    let parse_table = tuple((
        begin_bracket,
        table_line_names,
        table_line_id,
        line_ending,
        table_line_prot,
        table_line_s,
        newline,
        table_line_base1,
        newline,
        table_line_base2,
        newline,
        table_line_base3,
        space0,
        newline,
        space0,
        tag("}"),
    ));

    let parse_table = map(
        parse_table,
        |(
            _,
            vec_names,
            table_id,
            _,
            (_, _, _, _, prot, _, _, _, _),
            // _,
            (_, _, _, _starts),
            _,
            (_, _, _, base1),
            _,
            (_, _, _, base2),
            _,
            (_, _, _, base3),
            _,
            _,
            _,
            _,
        )| {
            // dbg!(prot);
            let shortname = if vec_names.len() > 1 {
                Some(vec_names[1].4.to_string())
            } else {
                None
            };

            let mut hash_codons = MultiMap::new();

            for idx in 0..prot.len() {
                let codon = vec![
                    base1.as_bytes()[idx],
                    base2.as_bytes()[idx],
                    base3.as_bytes()[idx],
                ];

                hash_codons.insert(prot.as_bytes()[idx], codon);
            }

            GeneticCodeTable {
                name: (vec_names[0].4).to_string(),
                shortname: shortname,
                id: table_id.3.parse::<u64>().unwrap(),
                start_codons: Vec::new(),
                end_codons: Vec::new(),
                table: hash_codons,
            }
        },
    );

    let (input, par) = separated_list1(
        terminated(tuple((space0, tag(","), space0)), opt(line_ending)),
        parse_table,
    )(input)?;

    // dbg!(par);

    Ok((input, par))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_open() {
        let fblastout = include_str!("../gc2.prt");

        let (_inp, res) = parse_input(fblastout).unwrap();

        for tb in res {
            let amino = tb.table.get_vec(&b'F').unwrap();
            for c in amino {
                let str_codon = std::str::from_utf8(c).unwrap();
                println!("{:?}", &str_codon);
            }
        }

        // dbg!(inp);

        // let res = Reader::from_file(path);
        assert_eq!(1, 1);
        // assert_eq!(res.unwrap().n_records, 2);
    }

    // #[test]
    // fn test_iter() {
    //     let fblastout = include_str!("../output_diamond.txt");

    //     let (inp, res) = parse_input(fblastout).unwrap();

    //     // let count_possible_frameshifts = res
    //     //     .into_iter()
    //     //     .tuple_windows()
    //     //     .filter(|(s1, s2)| s1.sseqid == s2.sseqid)
    //     //     .count();

    //     let count_possible_frameshifts = res
    //         .into_iter()
    //         .tuple_windows()
    //         .filter(|(s1, s2)| (s1.sseqid == s2.sseqid) & (s2.qframe == 2))
    //         .collect::<Vec<(BlastFeature, BlastFeature)>>();

    //     dbg!(count_possible_frameshifts);

    //     assert_eq!(1, 1);
    // }
}
