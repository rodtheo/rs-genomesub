// Inspired on https://pubmed.ncbi.nlm.nih.gov/34037690/

use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::alignment::AlignmentOperation::*;
use bio::data_structures::rank_select::RankSelect;
use bv::*;

use bio::data_structures::annot_map::AnnotMap;
use bio_types::annot::pos::Pos;
use bio_types::strand::ReqStrand;

use itertools::Itertools;
use itertools_num::ItertoolsNum;
// use bio::alignment::AlignmentOperation::*;
// use bio::alphabets::dna;
// use bio::io::fasta::{self, IndexedReader, Record};
// use bio::alignment::sparse::*;

#[derive(Debug, Clone)]
pub struct LiftBlock {
    // y is the reference or template sequence
    // x is the query
    pub alignment_op: Vec<AlignmentOperation>,
    // bitvectors to represent pairwise alignment between sequences y and x.
    // Where a 1 in bv_insertion indicates a insertion wrt reference
    // Where a 1 in bv_deletion indicates a deletion wrt reference
    pub bv_insertion: BitVec<u8>,
    pub bv_deletion: BitVec<u8>,
}

pub struct LiftBlockBuild {
    pub rs_bv2: RankSelect,
    pub rs_bv1: RankSelect,
}

impl LiftBlockBuild {
    pub fn new(liftblock: &LiftBlock) -> Self {
        LiftBlockBuild {
            rs_bv2: RankSelect::new(liftblock.bv_deletion.clone(), 1),
            rs_bv1: RankSelect::new(liftblock.bv_insertion.clone(), 1),
        }
    }

    pub fn liftover2assembly(&self, pos: u64) -> Option<u64> {
        // let bv2_rs = RankSelect::new(self.bv_deletion.clone(), 1);
        let bv2_rs = &self.rs_bv2;

        // let bv1_rs = RankSelect::new(self.bv_insertion.clone(), 1);
        let bv1_rs = &self.rs_bv1;

        // According to levioSAM, a position p in assembly matches a coordinate c in reference sequence
        // after we perform a select_0(p) with respect to (wrt) deletion-representation bitvector (bv2_rs).
        // Then this value is used as input to rank0 query wrt bv1_rs to obtain the corresponding position
        // on the reference sequence

        let sel = bv1_rs
            .select_0(pos + 1)
            .expect("Bitvector representing deletions is out of bounds.");

        let corresponding_pos = bv2_rs.rank_0(sel)?.checked_sub(1);

        corresponding_pos
    }
}

impl LiftBlock {
    pub fn new(pwalignment: &Vec<AlignmentOperation>) -> Self {
        let mut bv_in: BitVec<u8> = BitVec::new();
        let mut bv_del: BitVec<u8> = BitVec::new();

        for (i, op) in pwalignment.iter().enumerate() {
            if op == &Ins {
                bv_in.push(true);
            } else {
                bv_in.push(false);
            }

            if op == &Del {
                bv_del.push(true);
            } else {
                bv_del.push(false);
            }
        }

        LiftBlock {
            alignment_op: pwalignment.to_vec(),
            bv_insertion: bv_in,
            bv_deletion: bv_del,
        }
    }

    pub fn empty_matches(len: usize) -> Self {
        LiftBlock {
            alignment_op: vec![Match; len],
            bv_insertion: BitVec::new_fill(false, len as u64),
            bv_deletion: BitVec::new_fill(false, len as u64),
        }
    }
}

pub trait LiftOver {
    fn liftover2ref(&self, pos: u64) -> Option<u64>;
    fn liftover2ref_empower_rs(
        &self,
        pos: u64,
        bv_rs_deletion: &RankSelect,
        bv_rs_insertion: &RankSelect,
    ) -> Option<u64>;
    fn build_rankselects(&self) -> (RankSelect, RankSelect);
    fn liftover2assembly(&self, pos: u64) -> Option<u64>;
    fn liftover2assembly_empower_rs(
        &self,
        pos: u64,
        bv_rs_deletion: &RankSelect,
        bv_rs_insertion: &RankSelect,
    ) -> Option<u64>;
    fn concat_blocks(&self, other: LiftBlock) -> LiftBlock;
    fn insert_block(&self, other: LiftBlock, start_pos: u64, end_pos: u64) -> LiftBlock;
}

impl LiftOver for LiftBlock {
    // Note that pos is 0-based
    fn liftover2ref(&self, pos: u64) -> Option<u64> {
        let bv2_rs = RankSelect::new(self.bv_deletion.clone(), 1);
        let bv1_rs = RankSelect::new(self.bv_insertion.clone(), 1);

        // According to levioSAM, a position p in assembly matches a coordinate c in reference sequence
        // after we perform a select_0(p) with respect to (wrt) deletion-representation bitvector (bv2_rs).
        // Then this value is used as input to rank0 query wrt bv1_rs to obtain the corresponding position
        // on the reference sequence

        let sel = bv2_rs
            .select_0(pos + 1)
            .expect("Bitvector representing deletions is out of bounds.");

        let corresponding_pos = bv1_rs.rank_0(sel)?.checked_sub(1);

        corresponding_pos
    }

    fn liftover2assembly(&self, pos: u64) -> Option<u64> {
        let bv2_rs = RankSelect::new(self.bv_deletion.clone(), 1);
        let bv1_rs = RankSelect::new(self.bv_insertion.clone(), 1);

        // According to levioSAM, a position p in assembly matches a coordinate c in reference sequence
        // after we perform a select_0(p) with respect to (wrt) deletion-representation bitvector (bv2_rs).
        // Then this value is used as input to rank0 query wrt bv1_rs to obtain the corresponding position
        // on the reference sequence

        let sel = bv1_rs
            .select_0(pos + 1)
            .expect("Bitvector representing deletions is out of bounds.");

        let corresponding_pos = bv2_rs.rank_0(sel)?.checked_sub(1);

        corresponding_pos
    }

    fn concat_blocks(&self, other: LiftBlock) -> LiftBlock {
        let mut new_alignment_op = self.alignment_op.clone();
        new_alignment_op.extend(other.alignment_op);

        LiftBlock {
            alignment_op: new_alignment_op,

            bv_insertion: self
                .bv_insertion
                .bit_concat(other.bv_insertion)
                .to_bit_vec(),

            bv_deletion: self.bv_deletion.bit_concat(other.bv_deletion).to_bit_vec(),
        }
    }

    fn insert_block(&self, other: LiftBlock, start_pos: u64, end_pos: u64) -> LiftBlock {
        let big_block_insertion = self.bv_insertion.clone();

        // 5'prime bit vec
        let first_part_big_block = big_block_insertion.as_slice().bit_slice(..start_pos);

        // substitute pos 6,7,8 (0-based) with an array of trues
        let block_to_be_inserted = other.bv_insertion;

        // 3'prime bit vec
        let second_part_big_block = big_block_insertion.as_slice().bit_slice(end_pos..);

        let resulting_block_insertion = first_part_big_block
            .bit_concat(block_to_be_inserted.as_slice())
            .to_bit_vec()
            .bit_concat(second_part_big_block)
            .to_bit_vec();

        // doing the same for bv_deletion
        let big_block_deletion = self.bv_deletion.clone();

        // 5'prime bit vec - bv_deletion
        let first_part_big_block_deletion = big_block_deletion.as_slice().bit_slice(..start_pos);

        // substitute pos 6,7,8 (0-based) with an array of trues
        let block_to_be_inserted_deletion = other.bv_deletion;

        // 3'prime bit vec
        let second_part_big_block_deletion = big_block_deletion.as_slice().bit_slice(end_pos..);

        let resulting_block_deletion = first_part_big_block_deletion
            .bit_concat(block_to_be_inserted_deletion.as_slice())
            .to_bit_vec()
            .bit_concat(second_part_big_block_deletion)
            .to_bit_vec();

        let new_alignment_op = [
            &self.alignment_op[..start_pos as usize],
            &other.alignment_op,
            &self.alignment_op[end_pos as usize..],
        ]
        .concat();

        LiftBlock {
            alignment_op: new_alignment_op,
            bv_insertion: resulting_block_insertion,
            bv_deletion: resulting_block_deletion,
        }
    }

    fn liftover2ref_empower_rs(
        &self,
        pos: u64,
        bv_rs_deletion: &RankSelect,
        bv_rs_insertion: &RankSelect,
    ) -> Option<u64> {
        let bv2_rs = bv_rs_deletion;
        let bv1_rs = bv_rs_insertion;

        // According to levioSAM, a position p in assembly matches a coordinate c in reference sequence
        // after we perform a select_0(p) with respect to (wrt) deletion-representation bitvector (bv2_rs).
        // Then this value is used as input to rank0 query wrt bv1_rs to obtain the corresponding position
        // on the reference sequence

        let sel = bv2_rs
            .select_0(pos + 1)
            .expect("Bitvector representing deletions is out of bounds.");

        let corresponding_pos = bv1_rs.rank_0(sel)?.checked_sub(1);

        corresponding_pos
    }

    fn liftover2assembly_empower_rs(
        &self,
        pos: u64,
        bv_rs_deletion: &RankSelect,
        bv_rs_insertion: &RankSelect,
    ) -> Option<u64> {
        let bv2_rs = bv_rs_deletion;
        let bv1_rs = bv_rs_insertion;

        // According to levioSAM, a position p in assembly matches a coordinate c in reference sequence
        // after we perform a select_0(p) with respect to (wrt) deletion-representation bitvector (bv2_rs).
        // Then this value is used as input to rank0 query wrt bv1_rs to obtain the corresponding position
        // on the reference sequence

        let sel = bv1_rs
            .select_0(pos + 1)
            .expect("Bitvector representing deletions is out of bounds.");

        let corresponding_pos = bv2_rs.rank_0(sel)?.checked_sub(1);

        corresponding_pos
    }

    fn build_rankselects(&self) -> (RankSelect, RankSelect) {
        let bv2_rs = RankSelect::new(self.bv_deletion.clone(), 1);
        let bv1_rs = RankSelect::new(self.bv_insertion.clone(), 1);

        (bv2_rs, bv1_rs)
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, fs::File};

    use assert_cmd::assert::Assert;
    use itertools_num::ItertoolsNum;

    use crate::chain;

    use super::*;

    #[test]
    fn test_build_bv_alignment() {
        // let x = b"ACCGTGGAT";
        // let y = b"AAAAACCGTTGAT";
        let y = b"GCGATCGAGTGTACA";
        let x = b"GCTCCAGGCATGTACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block = LiftBlock::new(&alignment.operations);

        let p: u64 = 11;
        let lifted_pos = block.liftover2ref(p);

        assert_eq!(block.liftover2ref(0).unwrap(), 0);
        assert_eq!(block.liftover2ref(1).unwrap(), 1);
        assert_eq!(block.liftover2ref(2).unwrap(), 4);
        assert_eq!(block.liftover2ref(3).unwrap(), 5);
        assert_eq!(block.liftover2ref(4).unwrap(), 6);
        assert_eq!(block.liftover2ref(5).unwrap(), 7);
        assert_eq!(block.liftover2ref(6).unwrap(), 8);
        assert_eq!(block.liftover2ref(10).unwrap(), 9);
        assert_eq!(block.liftover2ref(11).unwrap(), 10);
        assert_eq!(block.liftover2ref(12).unwrap(), 11);
        assert_eq!(block.liftover2ref(15).unwrap(), 14);
    }

    #[test]
    fn test_build_bv_alignment_custom_rs() {
        // let x = b"ACCGTGGAT";
        // let y = b"AAAAACCGTTGAT";
        let y = b"GCGATCGAGTGTACA";
        let x = b"GCTCCAGGCATGTACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block = LiftBlock::new(&alignment.operations);

        let p: u64 = 11;
        let lifted_pos = block.liftover2ref(p);

        let (bv_rs_deletion, bv_rs_insertion) = block.build_rankselects();

        let res = vec![0, 1, 4, 5, 6, 7, 8];
        for i in (0..6).collect::<Vec<u64>>().iter() {
            assert_eq!(
                block
                    .liftover2ref_empower_rs(*i, &bv_rs_deletion, &bv_rs_insertion)
                    .unwrap(),
                res[*i as usize]
            );
        }
    }

    #[test]
    fn test_lift_reverse() {
        // let x = b"ACCGTGGAT";
        // let y = b"AAAAACCGTTGAT";
        let y = b"GCGATCGAGTGTACA";
        let x = b"GCTCCAGGCATGTACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block = LiftBlock::new(&alignment.operations);

        let p: u64 = 10;
        let lifted_pos = block.liftover2assembly(p);

        assert_eq!(lifted_pos.unwrap(), 11);
    }

    #[test]
    fn test_null_bitvector() {
        let y = b"GCGATCGAGTGTACA";
        let x = b"GCGATCGAGTGTACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block = LiftBlock::new(&alignment.operations);

        let p: u64 = 1;
        let lifted_pos = block.liftover2assembly(p);

        assert_eq!(lifted_pos.unwrap(), 1);
    }

    #[test]
    fn test_join_bvs() {
        let first: BitVec<u8> = BitVec::new_fill(false, 1);
        let middle: BitVec<u8> = BitVec::new_fill(true, 2);
        let end: BitVec<u8> = BitVec::new_fill(true, 1);

        let first_middle = first.bit_concat(middle);

        let first_middle_end = first_middle.bit_concat(end);
        // let first_and_last = [first.as_slice(), last.as_slice()].concat();

        assert_eq!(
            first_middle_end.to_bit_vec(),
            bit_vec![false, true, true, true]
        );
    }

    #[test]
    fn test_join_blocks() {
        let y = b"GCGA";
        let x = b"GCTCA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block_first = LiftBlock::new(&alignment.operations);

        let block_second = LiftBlock::new(&alignment.operations);

        dbg!(block_first.concat_blocks(block_second));
    }

    #[test]
    fn test_modify_inner_block() {
        let big_block: BitVec<u8> = BitVec::new_fill(false, 10);

        // 5'prime bit vec
        let first_part_big_block = big_block.as_slice().bit_slice(..6);

        // substitute pos 6,7,8 (0-based) with an array of trues
        let block_to_be_inserted: BitVec<u8> = BitVec::new_fill(true, 3);

        // 3'prime bit vec
        let second_part_big_block = big_block.as_slice().bit_slice(9..);

        let resulting_block = first_part_big_block
            .bit_concat(block_to_be_inserted.as_slice())
            .to_bit_vec()
            .bit_concat(second_part_big_block)
            .to_bit_vec();

        dbg!(resulting_block);
    }

    #[test]
    fn test_modify_inner_struct() {
        let y = b"GCGA";
        let x = b"GCTCA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        // gap open score: -5, gap extension score: -1
        let mut aligner = Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = aligner.global(x, y);

        let block_first = LiftBlock::new(&alignment.operations);

        dbg!(&block_first);
        let ref_seq_block = LiftBlock::empty_matches(10);

        let res = ref_seq_block.insert_block(block_first, 6, 9);

        dbg!(res);
    }

    #[test]
    fn test_naive_alignment_chain() {
        use super::*;
        use std::process::Command;
        use which::which;

        let has_bcftools = which("bcftools");

        match has_bcftools {
            Ok(_) => {
                let output_cmd_first = Command::new("tabix")
                    .arg("-p")
                    .arg("vcf")
                    .arg("examples/out.vcf.gz")
                    .output()
                    .expect("failed to execute tabix");

                let output_cmd = Command::new("bcftools")
                    .arg("consensus")
                    .arg("-f")
                    .arg("examples/input_diamond_sample_genome.fa")
                    .arg("-c")
                    .arg("examples/out.chain")
                    .arg("examples/out.vcf.gz")
                    .output()
                    .expect("failed to execute bcftools");

                dbg!(output_cmd);
            }
            Err(_) => {
                dbg!("Install bcftools and put it on executable PATH (http://samtools.github.io/bcftools/howtos/install.html)!");
                assert_eq!(0, 1)
            }
        }

        let chainout = include_str!("../out.chain");

        let (_, res) = chain::parse_file(chainout).unwrap();

        let chains: AnnotMap<String, LiftBlockBuild> = res.to_lift();

        let test_coords = vec![
            8976, 13677, 16288, 21533, 28095, 50907, 53599, 88945, 90288, 90282,
        ];
        let test_coords_res_pyliftover = vec![
            8975, 13677, 16290, 21534, 28098, 50911, 53605, 88952, 90296, 90281,
        ];

        for (idx, c) in test_coords.iter().enumerate() {
            let p: u64 = c - 1;
            let query = Pos::new("Pantoea_bnut".to_owned(), p as isize, ReqStrand::Forward);
            let hits = chains
                .find(&query)
                .map(|e| e.data())
                .collect::<Vec<&LiftBlockBuild>>();
            assert_eq!(hits.len(), 1);
            let lifted_pos = hits[0].liftover2assembly(p);
            dbg!(lifted_pos);
            // assert_eq!(lifted_pos.unwrap(), test_coords_res_pyliftover[idx]);
        }
    }

    #[test]
    fn test_lift_hg19_hg38() {
        // use super::*;
        use flate2::read::GzDecoder;
        use std::fs::File;
        use std::io::prelude::*;
        use std::path::Path;

        let p = Path::new("examples/hg19ToHg38.over.chain.gz");
        // dbg!(p);
        let mut reader = File::open(p).map(GzDecoder::new).unwrap();

        let mut data = String::new();
        reader.read_to_string(&mut data).unwrap();

        let (inp, res) = chain::parse_file(&data[..]).unwrap();

        let chains = res.to_lift();

        assert_eq!(res.n_records, 1278);

        let test_coords = vec![("chr8".to_string(), 141310715)];

        for (idx, (chr, c)) in test_coords.iter().enumerate() {
            let p: u64 = c - 1;
            let query = Pos::new(chr.clone(), p as isize, ReqStrand::Forward);
            let hits = chains
                .find(&query)
                .map(|e| e.data())
                .collect::<Vec<&LiftBlockBuild>>();
            assert_eq!(hits.len(), 1);
            let lifted_pos = hits[0].liftover2assembly(p);
            dbg!(lifted_pos);
            // assert_eq!(lifted_pos.unwrap(), test_coords_res_pyliftover[idx]);
        }
    }
}
