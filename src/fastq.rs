use std::io::{self, BufRead, BufReader, Read};

pub const PHRED_OFFSET: u8 = 33;
pub const MIN_PHRED: u8 = 2;

#[derive(Debug)]
pub struct FastqRecord {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

#[derive(Debug)]
pub enum FastqError {
    Io(io::Error),
    InvalidFormat(String),
    InvalidBase { pos: usize, byte: u8 },
    InvalidQuality { pos: usize, score: u8 },
    LengthMismatch { seq_len: usize, qual_len: usize },
}

impl From<io::Error> for FastqError {
    fn from(e: io::Error) -> Self {
        FastqError::Io(e)
    }
}

impl std::fmt::Display for FastqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FastqError::Io(e) => write!(f, "IO error: {}", e),
            FastqError::InvalidFormat(s) => write!(f, "invalid format: {}", s),
            FastqError::InvalidBase { pos, byte } => {
                write!(f, "invalid base at pos {}: 0x{:02x}", pos, byte)
            }
            FastqError::InvalidQuality { pos, score } => {
                write!(f, "invalid quality at pos {}: score={}", pos, score)
            }
            FastqError::LengthMismatch {
                seq_len,
                qual_len,
            } => write!(f, "length mismatch: seq={} qual={}", seq_len, qual_len),
        }
    }
}

#[inline]
pub fn encode_read_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(1),
        b'C' | b'c' => Some(2),
        b'G' | b'g' => Some(3),
        b'T' | b't' => Some(4),
        b'N' | b'n' => None,
        _ => None,
    }
}

#[inline]
pub fn validate_phred(qual_byte: u8) -> Option<u8> {
    let score = qual_byte.wrapping_sub(PHRED_OFFSET);
    if score >= MIN_PHRED && score <= 93 {
        Some(score)
    } else {
        None
    }
}

fn parse_record(lines: &[String]) -> Result<FastqRecord, FastqError> {
    if lines.len() != 4 {
        return Err(FastqError::InvalidFormat(format!(
            "expected 4 lines, got {}",
            lines.len()
        )));
    }
    if !lines[0].starts_with('@') {
        return Err(FastqError::InvalidFormat("header must start with @".into()));
    }
    if !lines[2].starts_with('+') {
        return Err(FastqError::InvalidFormat(
            "separator line must start with +".into(),
        ));
    }
    let name = lines[0][1..].split_whitespace().next().unwrap_or("").to_string();
    let seq_bytes = lines[1].as_bytes();
    let qual_bytes = lines[3].as_bytes();
    if seq_bytes.len() != qual_bytes.len() {
        return Err(FastqError::LengthMismatch {
            seq_len: seq_bytes.len(),
            qual_len: qual_bytes.len(),
        });
    }
    let mut seq = Vec::with_capacity(seq_bytes.len());
    let mut qual = Vec::with_capacity(qual_bytes.len());
    for (i, &b) in seq_bytes.iter().enumerate() {
        match encode_read_base(b) {
            Some(c) => seq.push(c),
            None => {
                return Err(FastqError::InvalidBase { pos: i, byte: b });
            }
        }
    }
    for (i, &q) in qual_bytes.iter().enumerate() {
        match validate_phred(q) {
            Some(s) => qual.push(s),
            None => {
                return Err(FastqError::InvalidQuality { pos: i, score: q });
            }
        }
    }
    Ok(FastqRecord { name, seq, qual })
}

pub struct FastqStreamReader<R: Read> {
    reader: BufReader<R>,
    line_buf: String,
    record_count: u64,
    error_count: u64,
}

impl<R: Read> FastqStreamReader<R> {
    pub fn new(reader: R) -> Self {
        FastqStreamReader {
            reader: BufReader::with_capacity(8 * 1024 * 1024, reader),
            line_buf: String::with_capacity(1024),
            record_count: 0,
            error_count: 0,
        }
    }

    pub fn stats(&self) -> (u64, u64) {
        (self.record_count, self.error_count)
    }

    fn read_line_trimmed(&mut self) -> io::Result<Option<String>> {
        self.line_buf.clear();
        let n = self.reader.read_line(&mut self.line_buf)?;
        if n == 0 {
            return Ok(None);
        }
        let trimmed = self.line_buf.trim().to_string();
        Ok(Some(trimmed))
    }

    pub fn next_record(&mut self) -> Option<Result<FastqRecord, FastqError>> {
        let mut lines: [String; 4] = [
            String::new(),
            String::new(),
            String::new(),
            String::new(),
        ];
        for i in 0..4 {
            match self.read_line_trimmed() {
                Ok(Some(line)) => {
                    if i == 0 && line.is_empty() {
                        return self.next_record();
                    }
                    lines[i] = line;
                }
                Ok(None) => {
                    if i == 0 {
                        return None;
                    }
                    self.error_count += 1;
                    return Some(Err(FastqError::InvalidFormat("truncated record".into())));
                }
                Err(e) => {
                    self.error_count += 1;
                    return Some(Err(FastqError::Io(e)));
                }
            }
        }
        self.record_count += 1;
        Some(parse_record(&lines))
    }
}

pub fn stream_fastq_file<P: AsRef<std::path::Path>>(
    path: P,
) -> io::Result<FastqStreamReader<std::fs::File>> {
    let file = std::fs::File::open(path)?;
    Ok(FastqStreamReader::new(file))
}
