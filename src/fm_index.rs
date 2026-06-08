use crate::bwt::BwtResult;
use crate::fasta::encode_base;
use crate::rank_select::{CheckpointedOcc, PackedBwt};
use serde::{Deserialize, Serialize};

const SA_SAMPLE_RATE: usize = 32;
const OCC_SAMPLE_RATE: usize = 128;

#[derive(Serialize, Deserialize)]
pub struct FmIndex {
    pub packed_bwt: PackedBwt,
    pub sa_samples: Vec<u32>,
    pub sa_sample_rate: usize,
    pub text_len: usize,
    pub seq_name: String,
    pub checkpointed_occ: CheckpointedOcc,
    pub sa_full: Vec<u32>,
}

impl FmIndex {
    pub fn build(text: &[u8], seq_name: &str) -> Self {
        let n = text.len();
        eprintln!("[fm-index] Building SA for {} bases...", n);
        let bwt_result = BwtResult::from_text(text);
        eprintln!("[fm-index] SA built. Building rank structures...");

        let packed_bwt = PackedBwt::new(&bwt_result.bwt, bwt_result.c_table);
        eprintln!("[fm-index] PackedBwt built. Building checkpointed OCC...");

        let checkpointed_occ = CheckpointedOcc::new(&bwt_result.bwt, bwt_result.c_table, OCC_SAMPLE_RATE);
        eprintln!("[fm-index] Checkpointed OCC built. Sampling SA...");

        let mut sa_samples = Vec::new();
        for i in (0..n).step_by(SA_SAMPLE_RATE) {
            sa_samples.push(bwt_result.sa[i]);
        }

        eprintln!(
            "[fm-index] Done. Text={} BWT={} SA_samples={}",
            n,
            packed_bwt.len(),
            sa_samples.len()
        );

        FmIndex {
            packed_bwt,
            sa_samples,
            sa_sample_rate: SA_SAMPLE_RATE,
            text_len: n,
            seq_name: seq_name.to_string(),
            checkpointed_occ,
            sa_full: bwt_result.sa,
        }
    }

    #[inline]
    pub fn lf_mapping(&self, i: usize) -> usize {
        self.checkpointed_occ.lf(i)
    }

    #[inline]
    pub fn lf_mapping_fast(&self, i: usize) -> usize {
        let c = self.checkpointed_occ.bwt[i];
        (self.checkpointed_occ.c_table[c as usize] + self.checkpointed_occ.occ(c, i)) as usize
    }

    pub fn backward_search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        let c_table = &self.checkpointed_occ.c_table;
        let occ = |c: u8, pos: usize| -> u32 { self.checkpointed_occ.occ(c, pos) };

        let m = pattern.len();
        if m == 0 {
            return Some((0, self.text_len));
        }
        let last = pattern[m - 1];
        let c = encode_base(last);
        if c == 0 {
            return None;
        }
        let mut lo = c_table[c as usize] as usize;
        let mut hi = c_table[c as usize] as usize;
        {
            let mut count = 0u32;
            for &b in &self.checkpointed_occ.bwt {
                if b == c {
                    count += 1;
                }
            }
            hi += count as usize;
        }
        if lo >= hi {
            return None;
        }

        for i in (0..m - 1).rev() {
            let ch = pattern[i];
            let c = encode_base(ch);
            if c == 0 {
                return None;
            }
            lo = c_table[c as usize] as usize + occ(c, lo) as usize;
            hi = c_table[c as usize] as usize + occ(c, hi) as usize;
            if lo >= hi {
                return None;
            }
        }
        Some((lo, hi))
    }

    pub fn backward_search_encoded(&self, encoded_pattern: &[u8]) -> Option<(usize, usize)> {
        let c_table = &self.checkpointed_occ.c_table;
        let occ = |c: u8, pos: usize| -> u32 { self.checkpointed_occ.occ(c, pos) };

        let m = encoded_pattern.len();
        if m == 0 {
            return Some((0, self.text_len));
        }
        let last = encoded_pattern[m - 1];
        if last == 0 {
            return None;
        }
        let mut lo = c_table[last as usize] as usize;
        let mut hi = c_table[last as usize] as usize;
        {
            let mut count = 0u32;
            for &b in &self.checkpointed_occ.bwt {
                if b == last {
                    count += 1;
                }
            }
            hi += count as usize;
        }
        if lo >= hi {
            return None;
        }

        for i in (0..m - 1).rev() {
            let c = encoded_pattern[i];
            if c == 0 {
                return None;
            }
            lo = c_table[c as usize] as usize + occ(c, lo) as usize;
            hi = c_table[c as usize] as usize + occ(c, hi) as usize;
            if lo >= hi {
                return None;
            }
        }
        Some((lo, hi))
    }

    pub fn resolve_offset(&self, sa_pos: usize) -> u32 {
        if sa_pos < self.sa_full.len() {
            return self.sa_full[sa_pos];
        }
        let mut pos = sa_pos;
        let mut steps = 0u32;
        loop {
            let sampled_idx = pos / self.sa_sample_rate;
            if pos % self.sa_sample_rate == 0 && sampled_idx < self.sa_samples.len() {
                let sa_val = self.sa_samples[sampled_idx];
                return sa_val.wrapping_add(steps) % self.text_len as u32;
            }
            pos = self.lf_mapping(pos);
            steps += 1;
            if steps as usize > self.text_len {
                break;
            }
        }
        u32::MAX
    }

    pub fn find_offsets(&self, pattern: &[u8]) -> Vec<u32> {
        let encoded: Vec<u8> = pattern.iter().map(|&b| encode_base(b)).collect();
        if encoded.iter().any(|&c| c == 0) {
            return vec![];
        }
        match self.backward_search_encoded(&encoded) {
            Some((lo, hi)) => {
                let mut offsets = Vec::with_capacity(hi - lo);
                for i in lo..hi {
                    offsets.push(self.resolve_offset(i));
                }
                offsets
            }
            None => vec![],
        }
    }

    pub fn find_offsets_encoded(&self, encoded_pattern: &[u8]) -> Vec<u32> {
        match self.backward_search_encoded(encoded_pattern) {
            Some((lo, hi)) => {
                let mut offsets = Vec::with_capacity(hi - lo);
                for i in lo..hi {
                    offsets.push(self.resolve_offset(i));
                }
                offsets
            }
            None => vec![],
        }
    }

    pub fn save<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        let writer = std::io::BufWriter::new(file);
        bincode::serialize_into(writer, self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }

    pub fn load<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        bincode::deserialize_from(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }
}
