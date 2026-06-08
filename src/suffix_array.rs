pub struct SuffixArray;

impl SuffixArray {
    pub fn build(text: &[u8]) -> Vec<u32> {
        let n = text.len();
        assert!(n > 0, "text must not be empty");
        assert!(n <= u32::MAX as usize, "text too large for u32 SA");
        if n <= 2 {
            return Self::build_naive(text);
        }
        Self::prefix_doubling(text)
    }

    fn build_naive(text: &[u8]) -> Vec<u32> {
        let n = text.len();
        let mut sa: Vec<u32> = (0..n as u32).collect();
        sa.sort_by(|&a, &b| {
            let ai = a as usize;
            let bi = b as usize;
            text[ai..].cmp(&text[bi..])
        });
        sa
    }

    fn prefix_doubling(text: &[u8]) -> Vec<u32> {
        let n = text.len();
        let mut sa: Vec<u32> = (0..n as u32).collect();
        let mut rank = vec![0u32; n];
        let mut tmp = vec![0u32; n];

        for i in 0..n {
            rank[i] = text[i] as u32;
        }

        let mut k = 1usize;
        while k < n {
            let rank_copy = rank.clone();
            let kk = k;

            sa.sort_by(|&a, &b| {
                let ai = a as usize;
                let bi = b as usize;
                let cmp0 = rank_copy[ai].cmp(&rank_copy[bi]);
                if cmp0 != std::cmp::Ordering::Equal {
                    return cmp0;
                }
                let a2 = if ai + kk < n { rank_copy[ai + kk] } else { 0 };
                let b2 = if bi + kk < n { rank_copy[bi + kk] } else { 0 };
                a2.cmp(&b2)
            });

            tmp[sa[0] as usize] = 1;
            for i in 1..n {
                let prev = sa[i - 1] as usize;
                let curr = sa[i] as usize;
                let same = rank_copy[prev] == rank_copy[curr]
                    && {
                        let rp = if prev + kk < n {
                            rank_copy[prev + kk]
                        } else {
                            0
                        };
                        let rc = if curr + kk < n {
                            rank_copy[curr + kk]
                        } else {
                            0
                        };
                        rp == rc
                    };
                tmp[curr] = tmp[prev] + if same { 0 } else { 1 };
            }

            std::mem::swap(&mut rank, &mut tmp);

            if rank[sa[n - 1] as usize] as usize == n {
                break;
            }
            k *= 2;
        }

        sa
    }
}
