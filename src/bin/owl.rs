use chrono::prelude::*;
use clap::{Args, Parser, Subcommand};
use indicatif::ProgressBar;
use log::debug;
use log::info;
use log::warn;
use log::LevelFilter;
use owl::all_rotations;
use owl::average_phred;
use owl::parse_custom_file_iter;
use owl::parse_repeat_file;
use owl::revcomp;
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

    /// Genomic region(s) to score for MSI (bed format with motif)
    #[arg(short, long)]
    regions: String,

    /// Sample name, SM bam tag is default.
    #[arg(short, long)]
    sample: Option<String>,

    /// Filter flag
    #[arg(short, long, default_value_t = 0)]
    flag: u16,

    /// Min MG
    #[arg(short, long, default_value_t = 0.98)]
    min_mg: f32,

    /// Max fraction of filtered reads
    #[arg(short, long, default_value_t = 0.2)]
    max_filt_frac: f64,

    /// Min coverage
    #[arg(short, long, default_value_t = 5)]
    min_cov: i32,

    /// Min mapQ
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

    /// Coefficient of variance cutoff
    #[arg(short, long, default_value_t = 5.0)]
    cov: f32,

    /// Exclude unphased reads
    #[arg(short, long)]
    unphased: bool,

    /// Canonicalize motifs (reverse complement, and rotate, keeping lowest)
    #[arg(short, long)]
    minimize: bool,

    /// Only score these samples
    #[arg(short, long, num_args = 1.., required = false)]
    samples: Vec<String>,

    /// Prefix for output files
    #[arg(short, long, required = true)]
    prefix: String,

    /// Min percentage sites for QC
    #[arg(short, long, default_value_t = 50.0)]
    min_pass: f32,

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
    ps: i32,
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
) -> (Region, Vec<ReadSegment>, f64, HashMap<i32, usize>) {
    let mut read_segs = Vec::new();
    let mut ps_sets: HashMap<i32, usize> = HashMap::new();

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
                "q-pos: {} {}, t-pos: {} {} mq: {} clip-l: {} clip-r: {}",
                r.pos(),
                r.reference_end(),
                t_start,
                t_end,
                r.mapq(),
                clips.0,
                clips.1,
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

        let ps = match r.aux(b"PS") {
            Ok(Aux::U8(val)) => val as i32,
            Ok(Aux::I8(val)) => val as i32,
            Ok(Aux::U16(val)) => val as i32,
            Ok(Aux::I16(val)) => val as i32,
            Ok(Aux::U32(val)) => val as i32,
            Ok(Aux::I32(val)) => val,
            _ => 0,
        };

        *ps_sets.entry(ps).or_insert(0) += 1;

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
            ps: ps,
        });
    }
    (rinfo, read_segs, nfail as f64 / ntotal as f64, ps_sets)
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
    "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"reference repeat length at this locus\">\n"
));
    header.push_str(&format!(
        "##INFO=<ID=MO,Number=1,Type=String,Description=\"repeat motif for this locus\">\n"
    ));
    header.push_str(&format!(
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"haplotype block\">\n"
    ));

    header.push_str(&format!(
        "##FORMAT=<ID=HP,Number=1,Type=Integer,Description=\"haplotype tag\">\n"
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

    let mut reader = IndexedReader::from_path(args.bam).expect("Failed to open BAM file");

    let mut ok_regions = 0;
    let mut bad_regions = 0;
    let mut region_count = 0;

    print!("{}", gen_header());
    println!("#Region\tinfo\tformat\t{}", samples.join(","));

    let pb = ProgressBar::new(regions.len() as u64);

    let mut total_phased = 0;
    let mut total_unphased = 0;

    for r in regions {
        region_count += 1;
        let control_reads = &get_reads(&r.0, &mut reader);

        pb.inc(1);

        let read_sub: (Region, Vec<ReadSegment>, f64, HashMap<i32, usize>) =
            process_region(control_reads, args.flag, &r.0, args.min_mg, args.min_mapq);

        let mut bad_region = false;

        if read_sub.2 > args.max_filt_frac {
            debug!(
                "Region {} failed because {:.2}% of reads failed filters.",
                r.0,
                read_sub.2 * 100.0,
            );
            bad_region = true;
            bad_regions += 1;
        };
        if control_reads.is_empty() {
            bad_regions += 1;
            debug!("Region {} contained no passing alignments.", r.0);
        }

        let (best_ps, _best_ps_count) = read_sub
            .3
            .iter()
            .filter(|(&ps, _)| ps != 0) // skip zero keys
            .max_by_key(|(_, &count)| count) // select max by count
            .map(|(&ps, &count)| (ps, count)) // extract pair
            .unwrap_or((0, 0)); // safe default

        let mut haps: HashMap<i32, Vec<(f64, String, String)>> = HashMap::new();

        for c in read_sub.1 {
            let aligner = MotifAligner::align(&c.seq, &r.1);

            let (score, start, end, pid) = aligner.backtrace();

            debug!(
                "hap:{} ps:{} seq-len:{} qual:{:.2} name:{} score/pid: {}/{} mlen:{} seq:{}",
                c.hap,
                c.ps,
                end - start, //
                average_phred(&c.qual),
                c.name,
                score,
                pid,
                r.1.len(),
                underline_span(&c.seq, start, end)
            );

            if average_phred(&c.qual) < args.min_bq
                || end - start < r.1.len()
                || pid < args.min_idt
                || (c.ps != 0 && c.ps != best_ps)
            {
                continue;
            }

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

            debug!("mu:{} sd:{} cv:{}", mu, sd, cov);

            let info_str = format!(
                "{},{},{},{:.1},{:.1},{}",
                best_ps,
                h.0,
                repeats.len(),
                mu,
                cov,
                jlen
            );

            if (h.1.len() as i32) >= args.min_cov && !bad_region {
                hap_info.push(HapInfo {
                    idx: h.0,
                    rest: info_str,
                });
            } else {
                hap_info.push(HapInfo {
                    idx: h.0,
                    rest: format!("{},{},{},.,.,.", best_ps, h.0, repeats.len()).to_string(),
                });
            }
        }

        if hap_info.is_empty() {
            hap_info.push(HapInfo {
                idx: 0,
                rest: "0,0,0,.,.,.".to_string(),
            });
        }

        hap_info.sort_by_key(|h| h.idx);
        let joined = hap_info
            .iter()
            .map(|h| h.to_string()) // requires Display for HapInfo
            .collect::<Vec<_>>()
            .join(";");
        println!(
            "{}\tRL={};MO={}\t{}\t{}",
            r.0,
            read_sub.0.end - read_sub.0.start,
            r.1,
            "PS:HP:CT:MU:CV:LN",
            joined
        );
        ok_regions += 1;
        if hap_info.len() > 1 {
            total_phased += 1;
        } else {
            total_unphased += 1;
        }
    }
    log::info!(
        "Passing region count: {}, Filtered region count: {}, Phased/unphased: {}/{}, Total region count: {}",
        ok_regions,
        bad_regions,
        total_phased,
        total_unphased,
        region_count,
    );

    if total_phased == 0 {
        log::warn!("No phased regions profiled, ensure the BAM file was phased.")
    }

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
            mem_kb as f64 / (1024.0 * 1024.0)
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

    let mut output_phased = args.prefix.clone();
    output_phased += ".owl-phased-sites.txt";

    let mut output_phased_summary = args.prefix.clone();
    output_phased_summary += ".owl-phased-summary.txt";

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

    let mut output_phased_fh = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_phased)
        .unwrap();

    let mut output_phased_summary_fh = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_phased_summary)
        .unwrap();

    // NOTE: extended header with per-sample passing stats
    output_scores_fh
        .write_all(
            b"#sample\t#high\t#low\t%high\t#phased\t%phased\t#sites\t#passing\t%passing\tqc\n",
        )
        .unwrap();

    output_motif_fh
        .write_all(b"motif\t#high\t#low\t%high\n")
        .unwrap();

    output_phased_fh
        .write_all(b"#sample\tregion\tphase_block\tmotif\thap1\thap2\n")
        .unwrap();
    output_phased_summary_fh
    .write_all(b"#sample:phase_block\tlow\th1-high\th2-high\tco-high\tsum\texpected_cohigh\tratio_obs_exp\n")
    .unwrap();

    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let start_time = Instant::now();
    let records = parse_custom_file_iter(&args.file).unwrap();

    // Per-sample high/low (existing)
    let mut sd_counts: HashMap<String, (usize, usize)> = HashMap::new(); // (hi, lo)

    // NEW: per-sample site totals and passing-site counts
    let mut sample_site_counts: HashMap<String, (usize, usize, usize)> = HashMap::new(); // (total_sites, passing_sites)

    let mut sample_cohigh_counts: HashMap<String, (usize, usize, usize, usize)> = HashMap::new(); // (h1+h2, h1, h2)

    // Per-motif counts (global)
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
            let mut site_passed = false; // did ANY allele pass at this site for this sample?
            let mut site_phased = false;

            // motif counts bucket for this site (global across samples)

            let mut motif_key = record.motif.clone();

            if args.minimize {
                let mut motifs = Vec::new();
                motifs.push(motif_key.clone());
                motifs.push(revcomp(&motif_key.clone()));
                let rotations = all_rotations(&motifs);
                motif_key = rotations.first().unwrap().clone();
            }

            let motifc = motif_counts.entry(motif_key).or_insert((0, 0));

            if alleles.len() > 1 {
                site_phased = true;
            }

            let mut high_hap = [0; 2];

            let mut n_phased_haps = 0;

            let mut phase_block = 0;

            for allele in alleles {
                // 1. If the haplotype has low read depth, skip
                // 2. If the haplotype is unphased (0) and skipping unphased
                // 3. If the haplotype is unphased, but there are other haplotypes (phased), skip

                phase_block = allele.phase_block;

                if allele.count < args.min_depth
                    || (allele.haplotype == 0 && args.unphased)
                    || alleles.len() > 1 && allele.haplotype == 0
                {
                    continue;
                }
                site_passed = true;

                // cov is CV Score, not coverage
                if allele.cov >= args.cov {
                    entry.0 += 1; // hi
                    motifc.0 += 1;
                } else {
                    entry.1 += 1; // lo
                    motifc.1 += 1;
                }

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
                if allele.haplotype != 0 {
                    n_phased_haps += 1;
                }

                if allele.haplotype == 0 || allele.cov < args.cov {
                    continue;
                }

                // this is a high score on either haplotype 1 or 2
                high_hap[(allele.haplotype - 1) as usize] = 1;
            }

            let new_name = format!("{}:{}", sample_name, phase_block);

            let co_high = sample_cohigh_counts
                .entry(new_name.clone())
                .or_insert((0, 0, 0, 0));

            if n_phased_haps == 2 {
                let v = high_hap[0] + high_hap[1];

                match v {
                    0 => {
                        co_high.0 += 1;
                    }
                    1 => {
                        if high_hap[0] == 1 {
                            co_high.1 += 1;
                        } else {
                            co_high.2 += 1;
                        }
                    }
                    2 => {
                        co_high.3 += 1;
                    }
                    _ => {
                        panic!("unexpected hap count v = {}", v);
                    }
                }
                writeln!(
                    output_phased_fh,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    sample_name, record.region, phase_block, record.motif, high_hap[0], high_hap[1]
                )
                .unwrap();
            }

            // Update per-sample site counters
            let site_entry = sample_site_counts
                .entry(sample_name.clone())
                .or_insert((0, 0, 0));
            site_entry.0 += 1; // saw this site for this sample
            if site_passed {
                site_entry.1 += 1; // at least one allele passed
            }
            if site_phased {
                site_entry.2 += 1;
            }
        }
    }

    // Stable output order by sample name
    let mut sorted_samples: Vec<_> = sd_counts.iter().collect();
    sorted_samples.sort_by_key(|(sample_name, _)| *sample_name);

    for (k, v) in sample_cohigh_counts {
        let low = v.0 as f64;
        let h1_high = v.1 as f64;
        let h2_high = v.2 as f64;
        let co_high = v.3 as f64;

        let n = low + h1_high + h2_high + co_high;

        // Expected co-high under independence
        let expected = if n > 0.0 {
            let p_x = (h1_high + co_high) / n;
            let p_y = (h2_high + co_high) / n;
            p_x * p_y * n
        } else {
            0.0
        };

        // observed/expected (0.0 if expected == 0)
        let ratio = if expected > 0.0 {
            co_high / expected
        } else {
            0.0
        };

        writeln!(
            output_phased_summary_fh,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}",
            k,
            v.0,
            v.1,
            v.2,
            v.3,
            (v.0 + v.1 + v.2 + v.3),
            expected,
            ratio,
        )
        .unwrap();
    }

    for (sample_name, (hi, lo)) in sorted_samples {
        let denom_hilo = *hi + *lo;
        let pct_high = if denom_hilo > 0 {
            100.0 * (*hi as f64 / denom_hilo as f64)
        } else {
            0.0
        };

        let (sites_total, sites_pass, sites_phased) = sample_site_counts
            .get(sample_name)
            .copied()
            .unwrap_or((0, 0, 0));

        // "If none are passing the value is zero"
        let pct_pass = if sites_total == 0 || sites_pass == 0 {
            0.0
        } else {
            100.0 * (sites_pass as f64 / sites_total as f64)
        };

        // "If none are passing the value is zero"
        let pct_phase = if sites_total == 0 || sites_phased == 0 {
            0.0
        } else {
            100.0 * (sites_phased as f64 / sites_total as f64)
        };

        let qc = if pct_pass < args.min_pass as f64 {
            "fail"
        } else {
            "pass"
        };

        writeln!(
            output_scores_fh,
            "{}\t{}\t{}\t{:.2}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}",
            sample_name,
            hi,
            lo,
            pct_high,
            sites_phased,
            pct_phase,
            sites_total,
            sites_pass,
            pct_pass,
            qc
        )
        .unwrap();
    }

    // Motif summary (global)
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
        let mem_kb = proc.memory(); // kilobytes
        info!(
            "Score completed in {:.2?}, memory used: {:.2} MB",
            duration,
            mem_kb as f64 / (1024.0 * 1024.0)
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
