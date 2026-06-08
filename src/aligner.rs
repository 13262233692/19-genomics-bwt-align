use crate::fastq::{FastqRecord, FastqStreamReader};
use crate::fm_index::FmIndex;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicU64, Ordering};

pub struct AlignResult {
    pub read_name: String,
    pub offsets: Vec<u32>,
}

pub struct AlignStats {
    pub total_reads: u64,
    pub matched_reads: u64,
    pub total_hits: u64,
    pub invalid_reads: u64,
}

pub struct BatchAligner<'a> {
    fm_index: &'a FmIndex,
}

impl<'a> BatchAligner<'a> {
    pub fn new(fm_index: &'a FmIndex) -> Self {
        BatchAligner { fm_index }
    }

    pub fn align_read(&self, record: &FastqRecord) -> AlignResult {
        let offsets = self.fm_index.find_offsets_encoded(&record.seq);
        AlignResult {
            read_name: record.name.clone(),
            offsets,
        }
    }

    pub fn align_reads_stream<R: std::io::Read, W: Write>(
        &self,
        reader: &mut FastqStreamReader<R>,
        writer: &mut BufWriter<W>,
    ) -> AlignStats {
        let total_reads = AtomicU64::new(0);
        let matched_reads = AtomicU64::new(0);
        let total_hits = AtomicU64::new(0);
        let invalid_reads = AtomicU64::new(0);

        let batch_size = 1024;

        loop {
            let mut batch: Vec<FastqRecord> = Vec::with_capacity(batch_size);
            for _ in 0..batch_size {
                match reader.next_record() {
                    Some(Ok(record)) => {
                        batch.push(record);
                    }
                    Some(Err(_)) => {
                        invalid_reads.fetch_add(1, Ordering::Relaxed);
                    }
                    None => break,
                }
            }
            if batch.is_empty() {
                break;
            }

            let results: Vec<AlignResult> = batch
                .par_iter()
                .map(|record| self.align_read(record))
                .collect();

            for result in &results {
                total_reads.fetch_add(1, Ordering::Relaxed);
                if !result.offsets.is_empty() {
                    matched_reads.fetch_add(1, Ordering::Relaxed);
                    total_hits.fetch_add(result.offsets.len() as u64, Ordering::Relaxed);
                    for &offset in &result.offsets {
                        let _ = writeln!(writer, "{}\t{}", result.read_name, offset);
                    }
                }
            }
        }

        AlignStats {
            total_reads: total_reads.load(Ordering::Relaxed),
            matched_reads: matched_reads.load(Ordering::Relaxed),
            total_hits: total_hits.load(Ordering::Relaxed),
            invalid_reads: invalid_reads.load(Ordering::Relaxed),
        }
    }
}

pub fn align_single_read(fm_index: &FmIndex, pattern: &[u8]) -> Vec<u32> {
    fm_index.find_offsets(pattern)
}
