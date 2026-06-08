use std::fs::File;
use std::io::{BufRead, BufReader, Result};
use std::path::Path;

pub const ALPHA: [u8; 5] = [b'$', b'A', b'C', b'G', b'T'];

#[inline]
pub fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 3,
        b'T' | b't' => 4,
        _ => 0,
    }
}

#[inline]
pub fn decode_base(c: u8) -> u8 {
    ALPHA[c as usize]
}

pub struct FastaRecord {
    pub name: String,
    pub seq: Vec<u8>,
}

pub struct FastaReader {
    reader: BufReader<File>,
}

impl FastaReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(FastaReader {
            reader: BufReader::with_capacity(64 * 1024, file),
        })
    }

    pub fn read_all(self) -> Result<Vec<FastaRecord>> {
        let mut records = Vec::new();
        let mut cur_name = String::new();
        let mut cur_seq: Vec<u8> = Vec::with_capacity(256 * 1024 * 1024);
        for line in self.reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('>') {
                if !cur_name.is_empty() {
                    records.push(FastaRecord {
                        name: cur_name.clone(),
                        seq: cur_seq.clone(),
                    });
                    cur_seq.clear();
                }
                cur_name = trimmed[1..].split_whitespace().next().unwrap_or("").to_string();
            } else {
                for b in trimmed.bytes() {
                    let c = encode_base(b);
                    if c > 0 {
                        cur_seq.push(c);
                    }
                }
            }
        }
        if !cur_name.is_empty() {
            records.push(FastaRecord {
                name: cur_name,
                seq: cur_seq,
            });
        }
        Ok(records)
    }
}

pub fn read_fasta_to_text<P: AsRef<Path>>(path: P) -> Result<(String, Vec<u8>)> {
    let records = FastaReader::open(path)?.read_all()?;
    if records.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "no records"));
    }
    let r = &records[0];
    let mut text = r.seq.clone();
    text.push(0);
    Ok((r.name.clone(), text))
}
