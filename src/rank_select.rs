use serde::{Deserialize, Serialize};

const SUPERBLOCK_BITS: usize = 512;
const BLOCK_BITS: usize = 64;

#[derive(Serialize, Deserialize)]
pub struct RankDict {
    bitmap: Vec<u64>,
    len: usize,
    superblock: Vec<u64>,
    block: Vec<u16>,
    popcounts: Vec<u8>,
}

impl RankDict {
    pub fn new() -> Self {
        let mut pc = vec![0u8; 256];
        for i in 0..256u16 {
            pc[i as usize] = i.count_ones() as u8;
        }
        RankDict {
            bitmap: Vec::new(),
            len: 0,
            superblock: Vec::new(),
            block: Vec::new(),
            popcounts: pc,
        }
    }

    pub fn from_bool_iter<I: Iterator<Item = bool>>(iter: I, len: usize) -> Self {
        let n_words = (len + 63) / 64;
        let mut bitmap = vec![0u64; n_words];
        for (i, b) in iter.enumerate() {
            if b {
                bitmap[i / 64] |= 1u64 << (i % 64);
            }
        }
        let mut pc = vec![0u8; 256];
        for i in 0..256u16 {
            pc[i as usize] = i.count_ones() as u8;
        }
        let mut rd = RankDict {
            bitmap,
            len,
            superblock: Vec::new(),
            block: Vec::new(),
            popcounts: pc,
        };
        rd.build_index();
        rd
    }

    fn build_index(&mut self) {
        let n_blocks = (self.len + BLOCK_BITS - 1) / BLOCK_BITS;
        let n_super = (self.len + SUPERBLOCK_BITS - 1) / SUPERBLOCK_BITS;
        self.superblock = vec![0u64; n_super + 1];
        self.block = vec![0u16; n_blocks + 1];

        let mut running = 0u64;
        for b in 0..n_blocks {
            if b % (SUPERBLOCK_BITS / BLOCK_BITS) == 0 {
                self.superblock[b / (SUPERBLOCK_BITS / BLOCK_BITS)] = running;
                self.block[b] = 0;
            } else {
                self.block[b] = (running - self.superblock[b / (SUPERBLOCK_BITS / BLOCK_BITS)]) as u16;
            }
            let start_word = b * (BLOCK_BITS / 64);
            let end_word = std::cmp::min(start_word + (BLOCK_BITS / 64), self.bitmap.len());
            for w in start_word..end_word {
                running += self.bitmap[w].count_ones() as u64;
            }
        }
        self.superblock[n_super] = running;
        self.block[n_blocks] = (running - self.superblock[n_super - 1]) as u16;
    }

    #[inline]
    pub fn rank1(&self, pos: usize) -> u64 {
        if pos == 0 {
            return 0;
        }
        let pos = std::cmp::min(pos, self.len);
        let super_idx = pos / SUPERBLOCK_BITS;
        let block_idx = pos / BLOCK_BITS;
        let mut r = self.superblock[super_idx] + self.block[block_idx] as u64;
        let bit_offset = block_idx * BLOCK_BITS;
        let word_offset = bit_offset / 64;
        let mut bits_left = pos - bit_offset;
        let mut w = word_offset;
        while bits_left >= 64 && w < self.bitmap.len() {
            r += self.bitmap[w].count_ones() as u64;
            bits_left -= 64;
            w += 1;
        }
        if bits_left > 0 && w < self.bitmap.len() {
            let mask = (1u64 << bits_left) - 1;
            r += (self.bitmap[w] & mask).count_ones() as u64;
        }
        r
    }

    #[inline]
    pub fn rank0(&self, pos: usize) -> u64 {
        pos as u64 - self.rank1(pos)
    }

    #[inline]
    pub fn access(&self, pos: usize) -> bool {
        if pos >= self.len {
            return false;
        }
        (self.bitmap[pos / 64] >> (pos % 64)) & 1 == 1
    }

    pub fn len(&self) -> usize {
        self.len
    }
}

#[derive(Serialize, Deserialize)]
pub struct WaveletTree {
    n: usize,
    sigma: u8,
    dicts: Vec<RankDict>,
    node_offsets: Vec<(usize, usize)>,
}

impl WaveletTree {
    pub fn new(data: &[u8], sigma: u8) -> Self {
        let n = data.len();
        if sigma <= 2 {
            let iter = data.iter().map(|&c| c == 1);
            let rd = RankDict::from_bool_iter(iter, n);
            return WaveletTree {
                n,
                sigma,
                dicts: vec![rd],
                node_offsets: vec![(0, sigma as usize)],
            };
        }
        let mut dicts = Vec::new();
        let mut node_offsets = Vec::new();
        Self::build_recursive(data, 0, sigma as usize, &mut dicts, &mut node_offsets);
        WaveletTree {
            n,
            sigma,
            dicts,
            node_offsets,
        }
    }

    fn build_recursive(
        data: &[u8],
        lo: usize,
        hi: usize,
        dicts: &mut Vec<RankDict>,
        node_offsets: &mut Vec<(usize, usize)>,
    ) -> usize {
        let node_idx = dicts.len();
        node_offsets.push((lo, hi));
        dicts.push(RankDict::new());

        if lo + 1 >= hi {
            return node_idx;
        }

        let mid = (lo + hi) / 2;
        let bitmap: Vec<bool> = data.iter().map(|&c| c as usize >= mid).collect();
        let rd = RankDict::from_bool_iter(bitmap.into_iter(), data.len());
        dicts[node_idx] = rd;

        let mut left_data = Vec::new();
        let mut right_data = Vec::new();
        for &c in data {
            if (c as usize) < mid {
                left_data.push(c);
            } else {
                right_data.push(c);
            }
        }

        if !left_data.is_empty() {
            Self::build_recursive(&left_data, lo, mid, dicts, node_offsets);
        }
        if !right_data.is_empty() {
            Self::build_recursive(&right_data, mid, hi, dicts, node_offsets);
        }

        node_idx
    }

    #[inline]
    pub fn access(&self, pos: usize) -> u8 {
        if self.sigma <= 2 {
            return if self.dicts[0].access(pos) { 1 } else { 0 };
        }
        let mut node = 0;
        let mut pos = pos;
        let mut lo = 0usize;
        let mut hi = self.sigma as usize;
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            let bit = self.dicts[node].access(pos);
            let r0 = self.dicts[node].rank0(pos + 1) as usize;
            let _r1 = self.dicts[node].rank1(pos + 1) as usize;
            if bit {
                pos = pos - r0;
                lo = mid;
                node = 2 * node + 2;
            } else {
                pos = r0 - 1;
                hi = mid;
                node = 2 * node + 1;
            }
            if node >= self.dicts.len() {
                break;
            }
        }
        lo as u8
    }

    #[inline]
    pub fn rank(&self, c: u8, pos: usize) -> u64 {
        if self.sigma <= 2 {
            if c == 0 {
                return self.dicts[0].rank0(pos);
            } else {
                return self.dicts[0].rank1(pos);
            }
        }
        let mut node = 0;
        let mut pos = pos;
        let mut lo = 0usize;
        let mut hi = self.sigma as usize;
        while lo + 1 < hi && node < self.dicts.len() {
            let mid = (lo + hi) / 2;
            if c as usize >= mid {
                let _r1 = self.dicts[node].rank1(pos);
                let r0 = self.dicts[node].rank0(pos);
                pos = (pos as u64 - r0) as usize;
                lo = mid;
                node = 2 * node + 2;
            } else {
                pos = self.dicts[node].rank0(pos) as usize;
                hi = mid;
                node = 2 * node + 1;
            }
        }
        pos as u64
    }
}

#[derive(Serialize, Deserialize)]
pub struct PackedBwt {
    pub n: usize,
    pub wavelet: WaveletTree,
    pub c_table: [u32; 5],
}

impl PackedBwt {
    pub fn new(bwt: &[u8], c_table: [u32; 5]) -> Self {
        let n = bwt.len();
        let sigma = 5u8;
        let wavelet = WaveletTree::new(bwt, sigma);
        PackedBwt {
            n,
            wavelet,
            c_table,
        }
    }

    #[inline]
    pub fn access(&self, pos: usize) -> u8 {
        self.wavelet.access(pos)
    }

    #[inline]
    pub fn rank(&self, c: u8, pos: usize) -> u64 {
        self.wavelet.rank(c, pos)
    }

    #[inline]
    pub fn lf(&self, i: usize) -> usize {
        let c = self.access(i);
        (self.c_table[c as usize] as u64 + self.rank(c, i)) as usize
    }

    pub fn len(&self) -> usize {
        self.n
    }
}

#[derive(Serialize, Deserialize)]
pub struct CheckpointedOcc {
    pub checkpoints: Vec<[u32; 5]>,
    pub sample_rate: usize,
    pub bwt: Vec<u8>,
    pub c_table: [u32; 5],
}

impl CheckpointedOcc {
    pub fn new(bwt: &[u8], c_table: [u32; 5], sample_rate: usize) -> Self {
        let n = bwt.len();
        let n_checkpoints = (n + sample_rate) / sample_rate + 1;
        let mut checkpoints = vec![[0u32; 5]; n_checkpoints];
        let mut counts = [0u32; 5];
        let mut cp_idx = 1usize;
        checkpoints[0] = [0; 5];
        for i in 0..n {
            counts[bwt[i] as usize] += 1;
            if (i + 1) % sample_rate == 0 && cp_idx < n_checkpoints {
                checkpoints[cp_idx] = counts;
                cp_idx += 1;
            }
        }
        if cp_idx < n_checkpoints {
            checkpoints[cp_idx] = counts;
        }
        CheckpointedOcc {
            checkpoints,
            sample_rate,
            bwt: bwt.to_vec(),
            c_table,
        }
    }

    #[inline]
    pub fn occ(&self, c: u8, pos: usize) -> u32 {
        if pos == 0 {
            return 0;
        }
        let pos = std::cmp::min(pos, self.bwt.len());
        let cp_idx = pos / self.sample_rate;
        let base = self.checkpoints[cp_idx][c as usize];
        let start = cp_idx * self.sample_rate;
        let mut count = base;
        for i in start..pos {
            if self.bwt[i] == c {
                count += 1;
            }
        }
        count
    }

    #[inline]
    pub fn lf(&self, i: usize) -> usize {
        let c = self.bwt[i];
        (self.c_table[c as usize] + self.occ(c, i)) as usize
    }
}
