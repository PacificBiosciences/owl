use chrono::prelude::*;
use clap::{Args, Parser, Subcommand};
use indicatif::ProgressBar;
use log::debug;
use log::info;
use log::warn;
use log::LevelFilter;
use owl::average_phred;
use owl::parse_custom_file_iter;
use owl::parse_repeat_file;
use owl::underline_span;
use owl::MotifAligner;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Record;
use rust_htslib::bam::{IndexedReader, Read, Reader};
use statrs::statistics::Statistics;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::time::Instant;
use sysinfo::System as SysInfoSystem;

#[derive(Parser)]
#[command(
    name = "owl",
    about = "Tools for microsatellite instability analysis for HiFi data",
    author,                // uses Cargo.toml authors
    version,               // uses Cargo.toml version
    long_about = None,
    // Show authors in --help:
    help_template = "\n{about}\n\n{usage-heading}\n{usage}\n\n{all-args}\n"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Profile a BAM file
    Profile(ProfileArgs),

    /// Merge multiple profiles
    Merge(MergeArgs),

    /// Score profiles
    Score(ScoreArgs),
}

#[derive(Args)]
struct ProfileArgs {
    /// Input BAM file
    #[arg(short, long)]
    bam: String,

    /// Genomic region to extract (e.g., chr:start-end)
    #[arg(short, long)]
    regions: String,

    /// Overwrite sample name, SM bam tag is used by default.
    #[arg(short, long)]
    sample: Option<String>,

    /// Filter flag
    #[arg(short, long, default_value_t = 0)]
    flag: u16,

    /// Min MG
    #[arg(short, long, default_value_t = 0.98)]
    min_mg: f32,

    /// Max Fraction of Filtered Reads
    #[arg(short, long, default_value_t = 0.2)]
    max_filt_frac: f64,

    /// Min Coverage
    #[arg(short, long, default_value_t = 5)]
    min_cov: i32,

    /// Min MapQ
    #[arg(short, long, default_value_t = 10)]
    min_mapq: u8,

    /// Min region avg base qual
    #[arg(short, long, default_value_t = 25.0)]
    min_bq: f64,

    /// Min motif alignment identity
    #[arg(short, long, default_value_t = 90.0)]
    min_idt: f64,

    /// Keep motif spans in the output in the ln field.
    #[arg(short, long, default_value_t = false)]
    keep_lengths: bool,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

#[derive(Args)]
struct MergeArgs {
    /// Input files to merge
    #[arg(short, long, num_args = 1.., required = true)]
    files: Vec<String>,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

#[derive(Args)]
struct ScoreArgs {
    /// Input file to score
    #[arg(short, long, required = true)]
    file: String,

    /// Min depth
    #[arg(short, long, default_value_t = 5)]
    min_depth: u32,

    /// Coefficient of Variance cutoff
    #[arg(short, long, default_value_t = 5.0)]
    cov: f32,

    /// Include unphased reads
    #[arg(short, long)]
    unphased: bool,

    /// Only score these samples
    #[arg(short, long, num_args = 1.., required = false)]
    samples: Vec<String>,

    /// Prefix for output files
    #[arg(short, long, required = true)]
    prefix: String,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

pub struct HapInfo {
    idx: i32,
    rest: String,
}

impl fmt::Display for HapInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.rest)
    }
}

pub struct ReadSegment {
    name: String,
    seq: String,
    qual: String,
    _len: usize,
    hap: i32,
}

pub struct Region {
    _name: String,
    start: i64,
    end: i64,
}

fn read_has_large_soft_clip(record: &Record, threshold: u32) -> (bool, bool) {
    let cigars = record.cigar();
    let start_soft_clip = cigars
        .first()
        .map(|c| matches!(c, Cigar::SoftClip(n) if *n > threshold))
        .unwrap_or(false);
    let end_soft_clip = cigars
        .last()
        .map(|c| matches!(c, Cigar::SoftClip(n) if *n > threshold))
        .unwrap_or(false);
    (start_soft_clip, end_soft_clip)
}

fn get_sample_names(bam_path: &str) -> Vec<String> {
    let bam = Reader::from_path(bam_path).expect("Failed to open BAM file");
    let header_text = String::from_utf8_lossy(bam.header().as_bytes());

    let mut samples = Vec::new();

    for line in header_text.lines() {
        if line.starts_with("@RG") {
            for field in line.split('\t') {
                if field.starts_with("SM:") {
                    let sample = field.trim_start_matches("SM:");
                    samples.push(sample.to_string());
                }
            }
        }
    }

    samples
}

pub fn compute_query_interval(
    cigar: &[Cigar],
    read_start: i64,
    region_start: i64,
    region_end: i64,
) -> (usize, usize) {
    let mut ref_pos = read_start;
    let mut query_pos = 0;

    let mut qstart = None;
    let mut qend = None;

    for c in cigar {
        let (cigar_len, consumes_query, consumes_ref) = match c {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => (*len as i64, true, true),
            Cigar::Ins(len) => (*len as i64, true, false),
            Cigar::Del(len) | Cigar::RefSkip(len) => (*len as i64, false, true),
            Cigar::SoftClip(len) | Cigar::Pad(len) => (*len as i64, true, false),
            Cigar::HardClip(_) => (0, false, false), // doesn't consume query or ref
        };

        if consumes_ref && ref_pos + cigar_len > region_start && qstart.is_none() {
            // start overlaps region_start
            let offset = (region_start - ref_pos).max(0);
            qstart = Some((query_pos + if consumes_query { offset } else { 0 }) as usize);
        }

        if consumes_ref && ref_pos < region_end && ref_pos + cigar_len >= region_end {
            // end overlaps region_end
            let offset = (region_end - ref_pos).max(0);
            qend = Some((query_pos + if consumes_query { offset } else { 0 }) as usize);
            break;
        }

        if consumes_ref {
            ref_pos += cigar_len;
        }
        if consumes_query {
            query_pos += cigar_len;
        }
    }

    (
        qstart.expect("region_start not found in CIGAR alignment"),
        qend.expect("region_end not found in CIGAR alignment"),
    )
}

fn get_reads(region: &str, reader: &mut IndexedReader) -> Vec<Record> {
    debug!("region {}", region);
    let mut records = Vec::new();

    reader.fetch(region).expect("Failed to fetch region");

    for r in reader.records() {
        match r {
            Ok(rec) => records.push(rec),
            Err(e) => {
                eprintln!("Error reading BAM record: {}", e);
                break;
            }
        }
    }

    records
}

pub fn process_region(
    reads: &Vec<Record>,
    flag: u16,
    region: &str,
    mg: f32,
    min_mq: u8,
) -> (Region, Vec<ReadSegment>, f64) {
    let mut read_segs = Vec::new();

    let parts: Vec<&str> = region.split(':').collect();
    if parts.len() != 2 {
        panic!("Invalid region format. Expected format: chr:start-end");
    }
    let range: Vec<&str> = parts[1].split('-').collect();
    if range.len() != 2 {
        panic!("Invalid region range. Expected format: start-end");
    }
    let t_name = parts[0];
    let mut t_start: i64 = range[0].parse().expect("Invalid start position");
    let mut t_end: i64 = range[1].parse().expect("Invalid end position");

    let rinfo = Region {
        _name: t_name.to_string(),
        start: t_start,
        end: t_end,
    };

    t_start -= 5;
    t_end += 5;

    let mut ntotal: i32 = 0;
    let mut nfail: i32 = 0;

    for r in reads {
        let name = String::from_utf8_lossy(r.qname()).to_string();

        // find reads with large clips (50bp or greater)
        let clips = read_has_large_soft_clip(r, 50);

        // Skip reads that don't span the region
        if r.pos() > t_start || r.reference_end() < t_end {
            continue;
        }
        ntotal += 1;

        let mg_ok = match r.aux(b"mg") {
            Ok(Aux::Float(val)) => {
                if val < mg {
                    debug!("Skipping read due to MG = {}", val);
                    false
                } else {
                    true
                }
            }
            _ => {
                warn!("Warning keeping read without MG");
                true
            }
        };

        if (r.flags() & flag) > 0 || r.mapq() < min_mq || clips.0 || clips.1 || !mg_ok {
            nfail += 1;
            debug!(
                "q-pos: {} {}, t-pos: {} {} mq: {}",
                r.pos(),
                r.reference_end(),
                t_start,
                t_end,
                r.mapq()
            );
            continue;
        }

        let hap = match r.aux(b"HP") {
            Ok(Aux::U8(val)) => val as i32,
            Ok(Aux::I8(val)) => val as i32,
            Ok(Aux::U16(val)) => val as i32,
            Ok(Aux::I16(val)) => val as i32,
            Ok(Aux::U32(val)) => val as i32,
            Ok(Aux::I32(val)) => val,
            _ => 0,
        };

        let (q_start, q_end) = compute_query_interval(&r.cigar(), r.pos() as i64, t_start, t_end);
        let seq = String::from_utf8_lossy(&r.seq().as_bytes()).to_string();
        let qual = r
            .qual()
            .iter()
            .map(|&q| (q + 33) as char)
            .collect::<String>();

        let dna_seq = seq[q_start as usize..q_end as usize].to_owned();
        let qual_seq = qual[q_start as usize..q_end as usize].to_owned();
        let dna_len = dna_seq.len();

        read_segs.push(ReadSegment {
            name: name,
            seq: dna_seq,
            qual: qual_seq,
            _len: dna_len,
            hap: hap,
        });
    }
    (rinfo, read_segs, nfail as f64 / ntotal as f64)
}

fn gen_header() -> String {
    let mut header = String::new();
    let today = Local::now().format("%Y%m%d").to_string();

    header.push_str(&format!(
        "##source=owl_v{}\n##fileDate={}\n",
        env!("CARGO_PKG_VERSION"),
        today
    ));
    header.push_str(&format!(
        "##FORMAT=<ID=HP,Number=1,Type=Integer,Description=\"haplotype tag\">\n"
    ));
    header.push_str(&format!(
        "##FORMAT=<ID=CT,Number=1,Type=Integer,Description=\"number of reads, post filters\">\n"
    ));
    header.push_str(&format!(
        "##FORMAT=<ID=CT,Number=1,Type=Integer,Description=\"number of reads, post filters\">\n"
    ));
    header.push_str(&format!(
        "##FORMAT=<ID=MU,Number=1,Type=Float,Description=\"mean repeat length\">\n"
    ));
    header.push_str(&format!("##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"coefficent of variance on repeat length\">\n"));

    header.push_str(&format!(
        "##FORMAT=<ID=LN,Number=.,Type=Int,Description=\"lengths of repeats\">\n"
    ));

    header
}

fn run_profile(args: ProfileArgs) {
    let start_time = Instant::now();

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let regions = parse_repeat_file(&args.regions).unwrap();

    let mut samples: Vec<String> = Vec::new();

    match &args.sample {
        Some(sample) => samples.push(sample.clone()),
        None => samples = get_sample_names(&args.bam),
    }

    let mut reader = IndexedReader::from_path(args.bam).unwrap();

    let mut ok_regions = 0;
    let mut bad_regions = 0;
    let mut region_count = 0;

    print!("{}", gen_header());
    println!("#Region\tlen\tmotif\tformat\t{}", samples.join(","));

    let pb = ProgressBar::new(regions.len() as u64);

    for r in regions {
        region_count += 1;
        let control_reads = &get_reads(&r.0, &mut reader);

        pb.inc(1);

        let read_sub = process_region(control_reads, args.flag, &r.0, args.min_mg, args.min_mapq);

        let mut bad_region = false;

        if read_sub.2 > args.max_filt_frac {
            debug!(
                "Region {} failed because {:.2}% of reads failed filters.",
                r.0,
                read_sub.2 * 100.0
            );
            bad_region = true;
            bad_regions += 1;
        };
        if control_reads.is_empty() {
            bad_regions += 1;
            debug!("Region {} contained no passing alignments.", r.0);
        }

        let mut haps: HashMap<i32, Vec<(f64, String, String)>> = HashMap::new();

        for c in read_sub.1 {
            let aligner = MotifAligner::align(&c.seq, &r.1);

            let (score, start, end, pid) = aligner.backtrace();

            if average_phred(&c.qual) < args.min_bq || end - start < r.1.len() || pid < args.min_idt
            {
                continue;
            }

            debug!(
                "hap:{} seq-len:{} qual:{:.2} name:{} score/pid: {}/{} mlen:{} seq:{}",
                c.hap,
                end - start, //
                average_phred(&c.qual),
                c.name,
                score,
                pid,
                r.1.len(),
                underline_span(&c.seq, start, end)
            );

            let read_repeat = c.seq[start..end].to_string();
            haps.entry(c.hap) // Replace `c.group` with your intended key
                .or_insert_with(Vec::new)
                .push(((end - start) as f64, read_repeat, c.name));
        }

        let mut hap_info = Vec::new();

        for mut h in haps {
            let mut repeats: Vec<String> = Vec::new();
            let mut repeat_lens: Vec<f64> = Vec::new();

            h.1.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            for a in &h.1 {
                repeats.push(a.1.clone());
                repeat_lens.push(a.0);
            }

            let mut jlen = ".".to_string();

            if args.keep_lengths {
                jlen = repeat_lens
                    .iter()
                    .map(|n| n.to_string())
                    .collect::<Vec<_>>()
                    .join(":");
            }

            let mu = (&repeat_lens).mean();
            let sd = (&repeat_lens).std_dev();
            let cov = (sd / mu) * 100.0;

            let info_str = format!("{},{},{:.1},{:.1},{}", h.0, repeats.len(), mu, cov, jlen);

            if (h.1.len() as i32) >= args.min_cov && !bad_region {
                hap_info.push(HapInfo {
                    idx: h.0,
                    rest: info_str,
                });
            } else {
                hap_info.push(HapInfo {
                    idx: h.0,
                    rest: format!("{},{},.,.,.", h.0, repeats.len()).to_string(),
                });
            }
        }

        if hap_info.is_empty() {
            hap_info.push(HapInfo {
                idx: 0,
                rest: "0,0,.,.,.".to_string(),
            });
        }

        hap_info.sort_by_key(|h| h.idx);
        let joined = hap_info
            .iter()
            .map(|h| h.to_string()) // requires Display for HapInfo
            .collect::<Vec<_>>()
            .join(";");
        println!(
            "{}\t{}\t{}\t{}\t{}",
            r.0,
            read_sub.0.end - read_sub.0.start,
            r.1,
            "hp:ct:mu:cv:ln",
            joined
        );
        ok_regions += 1;
    }
    log::info!(
        "Passing region count: {}, Filtered region count: {}, Total region count {}",
        ok_regions,
        bad_regions,
        region_count,
    );
    let duration = start_time.elapsed();
    let mut sys = SysInfoSystem::new_all();
    sys.refresh_processes();

    pb.finish_with_message("Hunted all regions!");

    if let Some(proc) = sys.process(sysinfo::get_current_pid().unwrap()) {
        let mem_kb = proc.memory(); // in kilobytes
        log::info!(
            "Profile completed in {:.2?}, memory used: {} MB",
            duration,
            mem_kb / (1024 * 1024)
        );
    } else {
        log::info!(
            "Profile completed in {:.2?}, memory usage unknown",
            duration
        );
    }
}

fn run_merge(args: MergeArgs) {
    let start_time = Instant::now();

    debug!("Merging {} files", args.files.len());

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    // Open all files and wrap them in BufReaders
    let mut readers: Vec<BufReader<File>> = args
        .files
        .iter()
        .map(|filename| {
            let file = File::open(Path::new(filename))
                .unwrap_or_else(|e| panic!("Failed to open {}: {}", filename, e));
            BufReader::new(file)
        })
        .collect();

    let mut eof = false;

    while !eof {
        let mut merged_line: Vec<String> = Vec::new();

        for (i, reader) in readers.iter_mut().enumerate() {
            let mut buf = String::new();
            let bytes = reader.read_line(&mut buf).unwrap();

            if bytes == 0 {
                eof = true;
                break;
            }

            let cols: Vec<&str> = buf.trim_end().split('\t').collect();
            if i == 0 {
                merged_line.extend(cols.iter().map(|s| s.to_string()));
            } else if cols.len() > 4 {
                merged_line.extend(cols[4..].iter().map(|s| s.to_string()));
            } else {
                // Pad with empty columns if needed
                //merged_line.extend(vec![""; 0].iter().map(|s| s.to_string()));
            }
        }

        if !eof {
            println!("{}", merged_line.join("\t"));
        }
    }

    let duration = start_time.elapsed();
    let mut sys = SysInfoSystem::new_all();
    sys.refresh_processes();

    if let Some(proc) = sys.process(sysinfo::get_current_pid().unwrap()) {
        let mem_kb = proc.memory(); // in kilobytes
        info!(
            "Merge completed in {:.2?}, memory used: {:.2} MB",
            duration,
            mem_kb as f64 / 1024.0
        );
    } else {
        info!("Merge completed in {:.2?}, memory usage unknown", duration);
    }
}

fn run_score(args: ScoreArgs) {
    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    let mut output_scores = args.prefix.clone();
    output_scores += ".owl-scores.txt";

    let mut output_motifs = args.prefix.clone();
    output_motifs += ".owl-motif-counts.txt";

    let mut output_scores_fh = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_scores)
        .unwrap();

    let mut output_motif_fh = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_motifs)
        .unwrap();

    output_scores_fh
        .write(b"#sample\t#high\t#low\t%high\n")
        .unwrap();
    output_motif_fh
        .write(b"motif\t#high\t#low\t%high\n")
        .unwrap();

    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let start_time = Instant::now();
    let records = parse_custom_file_iter(&args.file).unwrap();

    // Store running totals per sample
    let mut sd_counts: HashMap<String, (usize, usize)> = HashMap::new(); // (gt2, le2)
    let mut motif_counts: HashMap<String, (usize, usize)> = HashMap::new();

    let mut keepers = HashSet::new();
    for i in args.samples {
        keepers.insert(i);
    }

    for record_results in records {
        let record = record_results.unwrap();
        for (sample_name, alleles) in record.all_samples() {
            if !keepers.is_empty() && !keepers.contains(sample_name) {
                continue;
            }

            let entry = sd_counts.entry(sample_name.clone()).or_insert((0, 0));
            let motifc = motif_counts.entry(record.motif.clone()).or_insert((0, 0));
            for allele in alleles {
                if allele.count < args.min_depth || (allele.haplotype == 0 && !args.unphased) {
                    continue;
                }

                if allele.cov >= args.cov {
                    entry.0 += 1;
                    motifc.0 += 1;
                } else {
                    entry.1 += 1;
                    motifc.1 += 1;
                }

                // Optional: print per-allele info
                debug!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
                    record.region,
                    sample_name,
                    allele.haplotype,
                    allele.count,
                    allele.mean,
                    allele.cov,
                    allele.lengths
                );
            }
        }
    }

    let mut sorted_samples: Vec<_> = sd_counts.iter().collect();
    sorted_samples.sort_by_key(|(sample_name, _)| *sample_name);

    for (sample_name, (gt2, le2)) in sorted_samples {
        let frac = if *le2 > 0 {
            *gt2 as f64 / (*gt2 as f64 + *le2 as f64)
        } else {
            0.0
        };

        output_scores_fh
            .write(format!("{}\t{}\t{}\t{:.3}\n", sample_name, gt2, le2, frac * 100.0).as_bytes())
            .unwrap();
    }

    for (motif, (hi, lo)) in motif_counts {
        let total = hi + lo;
        let pct = if total == 0 {
            0.0
        } else {
            (hi as f64) / (total as f64) * 100.0
        };
        writeln!(output_motif_fh, "{motif}\t{hi}\t{lo}\t{pct:.2}").unwrap();
    }

    let duration = start_time.elapsed();
    let mut sys = SysInfoSystem::new_all();
    sys.refresh_processes();

    if let Some(proc) = sys.process(sysinfo::get_current_pid().unwrap()) {
        let mem_kb = proc.memory(); // in kilobytes
        info!(
            "Score completed in {:.2?}, memory used: {:.2} MB",
            duration,
            mem_kb as f64 / 1024.0
        );
    } else {
        info!("Score completed in {:.2?}, memory usage unknown", duration);
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Profile(args) => run_profile(args),
        Commands::Merge(args) => run_merge(args),
        Commands::Score(args) => run_score(args),
    }
}
