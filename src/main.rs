use clap::{Parser, Subcommand};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

mod aligner;
mod bwt;
mod fasta;
mod fastq;
mod fm_index;
mod rank_select;
mod suffix_array;

#[derive(Parser)]
#[command(name = "bwt-align", about = "FM-Index based genomic short-read aligner")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Build {
        #[arg(short, long, help = "Input FASTA file")]
        fasta: PathBuf,
        #[arg(short, long, help = "Output FM-Index file")]
        output: PathBuf,
    },
    Align {
        #[arg(short, long, help = "FM-Index file")]
        index: PathBuf,
        #[arg(short, long, help = "Input FASTQ file")]
        fastq: PathBuf,
        #[arg(short, long, help = "Output alignment file (default: stdout)")]
        output: Option<PathBuf>,
    },
    Search {
        #[arg(short, long, help = "FM-Index file")]
        index: PathBuf,
        #[arg(short, long, help = "Query sequence")]
        query: String,
    },
    Stats {
        #[arg(short, long, help = "FM-Index file")]
        index: PathBuf,
    },
}

fn cmd_build(fasta_path: &PathBuf, output_path: &PathBuf) {
    eprintln!("[build] Reading FASTA: {:?}", fasta_path);
    let t0 = Instant::now();
    let (seq_name, text) = fasta::read_fasta_to_text(fasta_path)
        .unwrap_or_else(|e| {
            eprintln!("[error] Failed to read FASTA: {}", e);
            std::process::exit(1);
        });
    eprintln!(
        "[build] Loaded {} bases ({:.1} MB) in {:.2}s",
        text.len(),
        text.len() as f64 / 1e6,
        t0.elapsed().as_secs_f64()
    );

    eprintln!("[build] Building FM-Index...");
    let t1 = Instant::now();
    let fm = fm_index::FmIndex::build(&text, &seq_name);
    eprintln!(
        "[build] FM-Index built in {:.2}s",
        t1.elapsed().as_secs_f64()
    );

    eprintln!("[build] Saving index to {:?}", output_path);
    let t2 = Instant::now();
    fm.save(output_path).unwrap_or_else(|e| {
        eprintln!("[error] Failed to save index: {}", e);
        std::process::exit(1);
    });
    let file_size = std::fs::metadata(output_path).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "[build] Index saved ({:.1} MB) in {:.2}s",
        file_size as f64 / 1e6,
        t2.elapsed().as_secs_f64()
    );
    eprintln!(
        "[build] Total time: {:.2}s",
        t0.elapsed().as_secs_f64()
    );
}

fn cmd_align(index_path: &PathBuf, fastq_path: &PathBuf, output_path: &Option<PathBuf>) {
    eprintln!("[align] Loading FM-Index from {:?}", index_path);
    let t0 = Instant::now();
    let fm = fm_index::FmIndex::load(index_path).unwrap_or_else(|e| {
        eprintln!("[error] Failed to load index: {}", e);
        std::process::exit(1);
    });
    eprintln!(
        "[align] Index loaded ({} bases) in {:.2}s",
        fm.text_len,
        t0.elapsed().as_secs_f64()
    );

    eprintln!("[align] Opening FASTQ: {:?}", fastq_path);
    let file = std::fs::File::open(fastq_path).unwrap_or_else(|e| {
        eprintln!("[error] Failed to open FASTQ: {}", e);
        std::process::exit(1);
    });
    let mut stream = fastq::FastqStreamReader::new(file);

    let out_writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(std::fs::File::create(p).unwrap_or_else(|e| {
            eprintln!("[error] Failed to create output: {}", e);
            std::process::exit(1);
        })),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(out_writer);

    let batch_aligner = aligner::BatchAligner::new(&fm);
    eprintln!("[align] Starting alignment with {} threads...", rayon::current_num_threads());
    let t1 = Instant::now();
    let stats = batch_aligner.align_reads_stream(&mut stream, &mut writer);
    let elapsed = t1.elapsed().as_secs_f64();

    eprintln!("[align] ────────────────────────────────────");
    eprintln!("[align] Alignment complete in {:.2}s", elapsed);
    eprintln!("[align] Total reads:    {}", stats.total_reads);
    eprintln!("[align] Matched reads:  {}", stats.matched_reads);
    eprintln!("[align] Total hits:     {}", stats.total_hits);
    eprintln!(
        "[align] Match rate:     {:.2}%",
        if stats.total_reads > 0 {
            stats.matched_reads as f64 / stats.total_reads as f64 * 100.0
        } else {
            0.0
        }
    );
    eprintln!(
        "[align] Throughput:     {:.0} reads/s",
        if elapsed > 0.0 {
            stats.total_reads as f64 / elapsed
        } else {
            0.0
        }
    );
    if stats.invalid_reads > 0 {
        eprintln!("[align] Invalid reads:  {}", stats.invalid_reads);
    }
    eprintln!("[align] ────────────────────────────────────");
}

fn cmd_search(index_path: &PathBuf, query: &str) {
    eprintln!("[search] Loading FM-Index from {:?}", index_path);
    let t0 = Instant::now();
    let fm = fm_index::FmIndex::load(index_path).unwrap_or_else(|e| {
        eprintln!("[error] Failed to load index: {}", e);
        std::process::exit(1);
    });
    eprintln!(
        "[search] Index loaded ({} bases) in {:.2}s",
        fm.text_len,
        t0.elapsed().as_secs_f64()
    );

    let t1 = Instant::now();
    let offsets = fm.find_offsets(query.as_bytes());
    let elapsed = t1.elapsed().as_secs_f64();

    if offsets.is_empty() {
        println!("No matches found for: {}", query);
    } else {
        println!("Found {} hit(s) for {}bp query in {:.4}s:", offsets.len(), query.len(), elapsed);
        for &off in &offsets {
            println!("  chr_offset: {}", off);
        }
    }
}

fn cmd_stats(index_path: &PathBuf) {
    eprintln!("[stats] Loading FM-Index from {:?}", index_path);
    let fm = fm_index::FmIndex::load(index_path).unwrap_or_else(|e| {
        eprintln!("[error] Failed to load index: {}", e);
        std::process::exit(1);
    });

    let file_size = std::fs::metadata(index_path).map(|m| m.len()).unwrap_or(0);
    println!("FM-Index Statistics");
    println!("─────────────────────────────");
    println!("Sequence name:    {}", fm.seq_name);
    println!("Text length:      {} bp", fm.text_len);
    println!("BWT length:       {}", fm.packed_bwt.len());
    println!("SA sample rate:   {}", fm.sa_sample_rate);
    println!("SA samples:       {}", fm.sa_samples.len());
    println!("SA full entries:  {}", fm.sa_full.len());
    println!("C-table:          {:?}", fm.checkpointed_occ.c_table);
    println!("OCC sample rate:  {}", fm.checkpointed_occ.sample_rate);
    println!("OCC checkpoints:  {}", fm.checkpointed_occ.checkpoints.len());
    println!("Index file size:  {:.1} MB", file_size as f64 / 1e6);
    println!(
        "Compression ratio: {:.2}x",
        fm.text_len as f64 / file_size as f64
    );
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Build { fasta, output } => cmd_build(&fasta, &output),
        Commands::Align {
            index,
            fastq,
            output,
        } => cmd_align(&index, &fastq, &output),
        Commands::Search { index, query } => cmd_search(&index, &query),
        Commands::Stats { index } => cmd_stats(&index),
    }
}
