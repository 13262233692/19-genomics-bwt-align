use crate::suffix_array::SuffixArray;

pub struct BwtResult {
    pub bwt: Vec<u8>,
    pub sa: Vec<u32>,
    pub c_table: [u32; 5],
}

impl BwtResult {
    pub fn from_text(text: &[u8]) -> Self {
        let n = text.len();
        let mut sa = vec![0u32; n];
        SuffixArray::build(text, &mut sa, n);
        let mut bwt = vec![0u8; n];
        for i in 0..n {
            let sai = sa[i] as usize;
            bwt[i] = if sai == 0 { text[n - 1] } else { text[sai - 1] };
        }

        let mut c_table = [0u32; 5];
        let mut freq = [0u32; 5];
        for &c in &bwt {
            freq[c as usize] += 1;
        }
        let mut acc = 0u32;
        for i in 0..5 {
            c_table[i] = acc;
            acc += freq[i];
        }

        BwtResult {
            bwt,
            sa,
            c_table,
        }
    }
}

pub fn compute_c_table(text: &[u8]) -> [u32; 5] {
    let mut freq = [0u32; 5];
    for &c in text {
        if (c as usize) < 5 {
            freq[c as usize] += 1;
        }
    }
    let mut c_table = [0u32; 5];
    let mut acc = 0u32;
    for i in 0..5 {
        c_table[i] = acc;
        acc += freq[i];
    }
    c_table
}
