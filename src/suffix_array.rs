pub struct SuffixArray;

const EMPTY: u32 = u32::MAX;

impl SuffixArray {
    pub fn build(text: &[u8], sa: &mut [u32], n: usize) {
        if n == 0 {
            return;
        }
        if n == 1 {
            sa[0] = 0;
            return;
        }
        let sigma = (*text.iter().max().unwrap() as usize).wrapping_add(1);
        if sigma > n || n <= 4 {
            let mut idx: Vec<u32> = (0..n as u32).collect();
            idx.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
            sa.copy_from_slice(&idx);
            return;
        }
        let mut bkt = vec![0u32; sigma];
        for &c in text.iter() {
            bkt[c as usize] += 1;
        }
        let mut sum = 0u32;
        for i in 0..sigma {
            let cnt = bkt[i];
            bkt[i] = sum;
            sum += cnt;
        }
        let bkt_head = bkt.clone();
        let mut bkt_tail = vec![0u32; sigma];
        let mut sum2 = 0u32;
        for i in 0..sigma {
            let old_head = bkt_head[i];
            let cnt = if i + 1 < sigma { bkt_head[i + 1] - old_head } else { n as u32 - old_head };
            sum2 += cnt;
            bkt_tail[i] = sum2;
        }

        let mut t = vec![false; n];
        t[n - 1] = true;
        for i in (0..n - 1).rev() {
            t[i] = if text[i] < text[i + 1] {
                true
            } else if text[i] > text[i + 1] {
                false
            } else {
                t[i + 1]
            };
        }

        for v in sa.iter_mut().take(n) {
            *v = EMPTY;
        }

        {
            let mut bkt = bkt_tail.clone();
            for i in (1..n).rev() {
                if t[i] && !t[i - 1] {
                    let c = text[i] as usize;
                    bkt[c] -= 1;
                    sa[bkt[c] as usize] = i as u32;
                }
            }
        }

        Self::induce_l(text, sa, &bkt_head, &t, n);
        Self::induce_s(text, sa, &bkt_tail, &t, n);

        let n1 = {
            let mut j = 0usize;
            for i in 0..n {
                let v = sa[i];
                if v != EMPTY && v > 0 && t[v as usize] && !t[v as usize - 1] {
                    sa[j] = v;
                    j += 1;
                }
            }
            for i in j..n {
                sa[i] = EMPTY;
            }
            j
        };

        if n1 == 0 {
            let mut idx: Vec<u32> = (0..n as u32).collect();
            idx.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
            sa.copy_from_slice(&idx);
            return;
        }

        let mut rank = vec![0u32; n];
        {
            let mut name = 0u32;
            let mut prev = sa[0] as usize;
            rank[prev] = 1;
            for i in 1..n1 {
                let curr = sa[i] as usize;
                let same = Self::lms_equal(text, &t, prev, curr, n);
                if !same {
                    name += 1;
                }
                rank[curr] = name + 1;
                prev = curr;
            }

            let all_unique = (name as usize + 1) == n1;

            let mut k = 0usize;
            for i in 0..n {
                if rank[i] > 0 {
                    sa[k] = rank[i] - 1;
                    k += 1;
                }
            }
            for i in k..n1 {
                sa[i] = 0;
            }

            if !all_unique {
                let max_name = name as usize + 1;
                if max_name <= 255 {
                    let mut sub_text = vec![0u8; n1 + 1];
                    for i in 0..n1 {
                        sub_text[i] = (sa[i] + 1) as u8;
                    }
                    sub_text[n1] = 0;
                    let mut sub_sa = vec![0u32; n1 + 1];
                    Self::build(&sub_text, &mut sub_sa, n1 + 1);
                    let mut j = 0usize;
                    for i in 0..=n1 {
                        if sub_sa[i] > 0 {
                            sa[j] = sub_sa[i] - 1;
                            j += 1;
                        }
                    }
                } else {
                    let mut sub_rank = vec![0u32; n1];
                    for i in 0..n1 {
                        sub_rank[i] = sa[i];
                    }
                    let mut sub_sa: Vec<u32> = (0..n1 as u32).collect();
                    let mut k = 1usize;
                    let mut tmp = vec![0u32; n1];
                    while k < n1 {
                        let kk = k;
                        sub_sa.sort_by(|&a, &b| {
                            let ra = sub_rank[a as usize];
                            let rb = sub_rank[b as usize];
                            match ra.cmp(&rb) {
                                std::cmp::Ordering::Equal => {
                                    let ra2 = if a as usize + kk < n1 { sub_rank[a as usize + kk] } else { 0 };
                                    let rb2 = if b as usize + kk < n1 { sub_rank[b as usize + kk] } else { 0 };
                                    ra2.cmp(&rb2)
                                }
                                ord => ord,
                            }
                        });
                        tmp[sub_sa[0] as usize] = 1;
                        for i in 1..n1 {
                            let prev = sub_sa[i - 1] as usize;
                            let curr = sub_sa[i] as usize;
                            let same = sub_rank[prev] == sub_rank[curr]
                                && {
                                    let rp = if prev + kk < n1 { sub_rank[prev + kk] } else { 0 };
                                    let rc = if curr + kk < n1 { sub_rank[curr + kk] } else { 0 };
                                    rp == rc
                                };
                            tmp[curr] = tmp[prev] + if same { 0 } else { 1 };
                        }
                        std::mem::swap(&mut sub_rank, &mut tmp);
                        if sub_rank[sub_sa[n1 - 1] as usize] as usize == n1 {
                            break;
                        }
                        k *= 2;
                    }
                    for i in 0..n1 {
                        sa[i] = sub_sa[i];
                    }
                }
            }
        }

        {
            let lms_pos: Vec<u32> = (1..n)
                .filter(|&i| t[i] && !t[i - 1])
                .map(|i| i as u32)
                .collect();
            let mut sorted_lms: Vec<u32> = vec![0; n1];
            for i in 0..n1 {
                let rank_val = sa[i] as usize;
                sorted_lms[i] = lms_pos[rank_val];
            }
            for v in sa.iter_mut().take(n) {
                *v = EMPTY;
            }
            {
                let mut bkt = bkt_tail.clone();
                for i in (0..n1).rev() {
                    let idx = sorted_lms[i] as usize;
                    let c = text[idx] as usize;
                    bkt[c] -= 1;
                    sa[bkt[c] as usize] = idx as u32;
                }
            }
        }

        Self::induce_l(text, sa, &bkt_head, &t, n);
        Self::induce_s(text, sa, &bkt_tail, &t, n);
    }

    fn lms_equal(text: &[u8], t: &[bool], a: usize, b: usize, n: usize) -> bool {
        if text[a] != text[b] {
            return false;
        }
        let mut i = 1usize;
        while a + i < n && b + i < n {
            if text[a + i] != text[b + i] {
                return false;
            }
            if i > 0 && ((t[a + i] && !t[a + i - 1]) || (t[b + i] && !t[b + i - 1])) {
                return (t[a + i] && !t[a + i - 1]) && (t[b + i] && !t[b + i - 1]);
            }
            i += 1;
        }
        false
    }

    fn induce_l(text: &[u8], sa: &mut [u32], bkt_head: &[u32], t: &[bool], n: usize) {
        let mut bkt = bkt_head.to_vec();
        for i in 0..n {
            if sa[i] == EMPTY || sa[i] == 0 {
                continue;
            }
            let j = (sa[i] - 1) as usize;
            if !t[j] {
                let c = text[j] as usize;
                sa[bkt[c] as usize] = j as u32;
                bkt[c] += 1;
            }
        }
    }

    fn induce_s(text: &[u8], sa: &mut [u32], bkt_tail: &[u32], t: &[bool], n: usize) {
        let mut bkt = bkt_tail.to_vec();
        for i in (0..n).rev() {
            if sa[i] == EMPTY || sa[i] == 0 {
                continue;
            }
            let j = (sa[i] - 1) as usize;
            if t[j] {
                let c = text[j] as usize;
                bkt[c] -= 1;
                sa[bkt[c] as usize] = j as u32;
            }
        }
    }
}
