use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Underlines a span of the input string using ANSI escape codes.
/// Returns a new `String` with the underlined range.
///
/// # Arguments
/// * `input` - The full string
/// * `start` - Start index of the span to underline (inclusive)
/// * `end` - End index of the span to underline (exclusive)
pub fn underline_span(input: &str, start: usize, end: usize) -> String {
    if start >= end || end > input.len() {
        return input.to_string(); // return unchanged if invalid range
    }

    let prefix = &input[..start];
    let span = &input[start..end];
    let suffix = &input[end..];

    format!("{}{}{}{}", prefix, "\x1b[4m", span, "\x1b[0m") + suffix
}

fn clean_motif(raw: &str) -> String {
    raw.trim_matches(|c| c == '(' || c == ')' || c == 'n')
        .to_string()
}

/// Parses the repeat file and returns a Vec of (region_string, motif)
pub fn parse_repeat_file(path: &str) -> io::Result<Vec<(String, String)>> {
    let reader = BufReader::new(File::open(path)?);
    let mut results = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() != 4 {
            eprintln!("Skipping malformed line: {line}");
            continue;
        }

        let chrom = fields[0];
        let start = fields[1];
        let end = fields[2];
        let motif = clean_motif(fields[3]);

        let region = format!("{chrom}:{start}-{end}");
        results.push((region, motif));
    }

    Ok(results)
}

/// Computes the average Phred quality score from a PacBio HiFi quality string.
/// Assumes Phred+33 ASCII encoding.
pub fn average_phred(phred_str: &str) -> f64 {
    if phred_str.is_empty() {
        return 0.0;
    }

    let total: u32 = phred_str
        .bytes()
        .map(|b| (b - 33) as u32) // Phred+33 decoding
        .sum();

    total as f64 / phred_str.len() as f64
}

#[derive(Debug, Clone, Copy)]
pub struct Cell {
    score: i32,
    i: usize,
    j: usize,
}

pub struct MotifAligner {
    matrix: Vec<Vec<Cell>>,
    target: Vec<u8>,
    motif: Vec<u8>,
}

impl MotifAligner {
    pub fn new(rows: usize, cols: usize) -> Self {
        let mut matrix = vec![
            vec![
                Cell {
                    score: 0,
                    i: 0,
                    j: 0,
                };
                cols
            ];
            rows
        ];

        // Initialize first column with gap penalties
        for i in 1..rows {
            matrix[i][0].score = 0;
            matrix[i][0].i = i - 1;
            matrix[i][0].j = 0;
        }

        // Initialize first row with gap penalties
        for j in 1..cols {
            matrix[0][j].score = 0;
            matrix[0][j].i = 0;
            matrix[0][j].j = j - 1;
        }

        Self {
            matrix,
            target: Vec::new(),
            motif: Vec::new(),
        }
    }

    pub fn set_cell(&mut self, i: usize, j: usize, score: i32, from_i: usize, from_j: usize) {
        self.matrix[i][j] = Cell {
            score,
            i: from_i,
            j: from_j,
        };
    }

    pub fn get_cell(&self, i: usize, j: usize) -> Cell {
        self.matrix[i][j]
    }

    pub fn align(target: &str, motif: &str) -> Self {
        let target_bytes = target.as_bytes().to_vec();
        let motif_bytes = motif.as_bytes().to_vec();
        let rows = target_bytes.len() + 1;
        let cols = motif_bytes.len() + 1;

        let mut aligner = Self::new(rows, cols);
        aligner.target = motif_bytes;
        aligner.motif = target_bytes;

        aligner.smith_waterman_align();

        aligner
    }

    fn smith_waterman_align(&mut self) {
        for i in 1..self.rows() {
            // First pass: compute all cells in the row
            for j in 1..self.cols() {
                let (score, from_i, from_j) = self.score_cell(i, j);
                self.set_cell(i, j, score, from_i, from_j);
            }

            // Check if last column score is greater than first column score
            let last = self.get_cell(i, self.cols() - 1);
            let first = self.get_cell(i, 0);

            if last.score > first.score {
                self.set_cell(i, 0, last.score, i, self.cols() - 1);
            }

            // Second pass: recompute starting from column 1 (not 0)
            for j in 1..self.cols() {
                let (score, from_i, from_j) = self.score_cell(i, j);
                self.set_cell(i, j, score, from_i, from_j);
            }
        }
    }

    fn score_cell(&self, i: usize, j: usize) -> (i32, usize, usize) {
        let match_score = 2;
        let mismatch_penalty = -2;
        let gap_penalty = -3;

        let diag = self.matrix[i - 1][j - 1].score
            + if self.motif[i - 1] == self.target[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };

        let up = self.matrix[i - 1][j].score + gap_penalty;
        let left = self.matrix[i][j - 1].score + gap_penalty;

        let (score, from_i, from_j) = {
            let mut best_score = 0;
            let mut src = (i, j); // default to self (zero)

            if diag > best_score {
                best_score = diag;
                src = (i - 1, j - 1);
            }
            if up > best_score {
                best_score = up;
                src = (i - 1, j);
            }
            if left > best_score {
                best_score = left;
                src = (i, j - 1);
            }

            (best_score, src.0, src.1)
        };

        (score, from_i, from_j)
    }

    pub fn best_score(&self) -> (i32, usize, usize) {
        let mut best = (0, 0, 0);
        for i in 0..self.rows() {
            for j in 0..self.cols() {
                let cell = self.matrix[i][j];
                if cell.score >= best.0 {
                    best = (cell.score, i, j);
                }
            }
        }
        best
    }

    pub fn rows(&self) -> usize {
        self.matrix.len()
    }

    pub fn cols(&self) -> usize {
        self.matrix[0].len()
    }

    pub fn print_matrix(&self) {
        // Print column indices
        print!("            ");
        for j in 0..self.cols() {
            print!("{:^18}", j);
        }
        println!();

        // Print column headers (sequence characters or '-')
        print!("            ");
        for j in 0..self.cols() {
            if j == 0 {
                print!("{:^18}", "-");
            } else {
                print!("{:^18}", self.target[j - 1] as char);
            }
        }
        println!();

        // Print each row with its index and character
        for i in 0..self.rows() {
            if i == 0 {
                print!("{:>2}  {:>2}  ", i, "-");
            } else {
                print!("{:>2}  {:>2}  ", i, self.motif[i - 1] as char);
            }

            for j in 0..self.cols() {
                let Cell {
                    score,
                    i: from_i,
                    j: from_j,
                } = self.matrix[i][j];
                print!("(i:{:>3} j:{:>3} s:{:>3}) ", from_i, from_j, score);
            }
            println!();
        }
    }
    /// Returns (score, motif_start, motif_end, percent_identity) from the best alignment.
    /// %identity = matches / (matches + mismatches + gaps) × 100
    pub fn backtrace(&self) -> (i32, usize, usize, f64) {
        let (score, mut i, mut j) = self.best_score();
        let end = i;

        let mut aligned_len: usize = 0; // counts all alignment columns (diag + gaps)
        let mut matches: usize = 0; // counts only diagonal matches

        loop {
            let cell = self.get_cell(i, j);
            // stop at local alignment boundary or self-pointer
            if (cell.i == i && cell.j == j) || cell.score == 0 {
                break;
            }

            let next_i = cell.i;
            let next_j = cell.j;
            let cols = self.cols();

            // Special "wrap-left" hop inserted by the DP (i,0) -> (i, cols-1):
            // this is bookkeeping, not a real alignment column, so don't count it.
            let is_wrap_left = i == next_i && j == 0 && next_j == cols - 1;

            if !is_wrap_left {
                // Diagonal: aligned pair -> contributes 1 column; may be match or mismatch
                if next_i + 1 == i && next_j + 1 == j {
                    aligned_len += 1;
                    if self.motif[i - 1] == self.target[j - 1] {
                        matches += 1;
                    }
                // Up: gap in motif (insertion in sequence) -> contributes 1 gap column
                } else if next_i + 1 == i && next_j == j {
                    aligned_len += 1;
                // Left: gap in sequence (insertion in motif) -> contributes 1 gap column
                } else if next_i == i && next_j + 1 == j {
                    aligned_len += 1;
                }
                // Any other pointer shape is ignored; with current DP this shouldn't occur.
            }

            i = next_i;
            j = next_j;
        }

        let start = i;
        let pid = if aligned_len > 0 {
            (matches as f64 / aligned_len as f64) * 100.0
        } else {
            0.0
        };

        (score, start, end, pid)
    }
}

#[derive(Debug)]
pub struct SampleData {
    pub haplotype: u32,
    pub count: u32,
    pub mean: f32,
    pub cov: f32,
    pub lengths: Vec<u32>,
}

#[derive(Debug)]
pub struct RegionRecord {
    pub region: String,
    pub len: u32,
    pub motif: String,
    pub sample_data: HashMap<String, Vec<SampleData>>,
}

impl RegionRecord {
    pub fn from_line(line: &str, sample_names: &[String]) -> Option<Self> {
        if line.trim().is_empty() || line.starts_with('#') {
            return None;
        }

        let mut cols = line.split_whitespace();
        let region = cols.next()?.to_string();
        let len: u32 = cols.next()?.parse().ok()?;
        let motif = cols.next()?.to_string();
        let format = cols.next()?.to_string();

        let remaining: Vec<&str> = cols.collect();
        /*
        if format_and_samples.len() < sample_names.len() * 2 {
            return None;
        }
        */

        let mut sample_data = HashMap::new();
        for (i, sample_name) in sample_names.iter().enumerate() {
            let data = remaining[i];
            let alleles = data
                .split(';')
                .filter_map(|s| parse_sample_entry(&format, s))
                .collect::<Vec<_>>();
            sample_data.insert(sample_name.clone(), alleles);
        }

        Some(Self {
            region,
            len,
            motif,
            sample_data,
        })
    }

    pub fn get_sample(&self, sample_name: &str) -> Option<&Vec<SampleData>> {
        self.sample_data.get(sample_name)
    }

    pub fn all_samples(&self) -> impl Iterator<Item = (&String, &Vec<SampleData>)> {
        self.sample_data.iter()
    }
}

fn parse_sample_entry(_format: &str, value: &str) -> Option<SampleData> {
    let fields: Vec<&str> = value.splitn(5, ',').collect();
    if fields.len() < 5 {
        return None;
    }

    let hp: u32 = fields[0].parse().ok()?;
    let ct: u32 = fields[1].parse().ok()?;
    let mu: f32 = fields[2].parse().ok()?;
    let cov: f32 = fields[3].parse().ok()?;
    let lengths: Vec<u32> = fields[4]
        .split(':')
        .filter_map(|s| s.parse().ok())
        .collect();

    Some(SampleData {
        haplotype: hp,
        count: ct,
        mean: mu,
        cov,
        lengths,
    })
}

pub struct RegionRecordIter<R: BufRead> {
    reader: R,
    header_fields: Vec<String>,
    line: String,
}

impl<R: BufRead> Iterator for RegionRecordIter<R> {
    type Item = io::Result<RegionRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line.clear();
        match self.reader.read_line(&mut self.line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                let line = self.line.trim_end().to_string();
                Some(
                    RegionRecord::from_line(&line, &self.header_fields)
                        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Parse error")),
                )
            }
            Err(e) => Some(Err(e)),
        }
    }
}

pub fn parse_custom_file_iter(path: &str) -> io::Result<RegionRecordIter<BufReader<File>>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut buf = String::new();

    // Produce header_fields exactly once; no unused assignment warning.
    let header_fields: Vec<String> = loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "No column header line starting with '#' found",
            ));
        }

        if buf.starts_with("##") {
            continue;
        }

        if buf.starts_with('#') {
            let cols_raw = buf.trim_start_matches('#').trim();

            let cols: Vec<&str> = if cols_raw.contains('\t') {
                cols_raw.split('\t').collect()
            } else {
                cols_raw.split_whitespace().collect()
            };

            if cols.len() < 5 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Malformed header line: {}", cols_raw),
                ));
            }

            let fields = cols.iter().skip(4).map(|s| s.to_string()).collect();
            break fields; // <- assign once via loop break
        }

        // else: non-header/non-metadata line before header; keep scanning
    };

    Ok(RegionRecordIter {
        reader,
        header_fields,
        line: String::new(),
    })
}

#[cfg(test)]

mod tests {
    use super::*;
    use assert_float_eq::assert_float_relative_eq;

    #[test]
    fn test_motif_alignment_basic() {
        let target = "ACGTTTTTTTTTTTTACGTTTTACGTACGTACGACGACGACGTTTTT";
        let motif = "ACG";
        let aligner = MotifAligner::align(target, motif);

        aligner.print_matrix();

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "motif:{} score:{} tstart:{} tend:{} pid:{:.3} tannoated:{}",
            motif,
            bt.0,
            bt.1,
            bt.2,
            bt.3,
            underline_span(target, bt.1, bt.2)
        );

        assert_eq!(score, 30);
        assert_eq!(bt.1, 22);
        assert_eq!(bt.2, 42);
        assert_eq!(bt.3, 90.0);
    }

    #[test]
    fn test_motif_alignment_real() {
        let target =
            "TCTTCTTCTTTTCTTTTTTTTCTTTTCTTTCCTTTTCCTTTTCTTTTCTTTTCTTTTCTTTTCTCTTCTCTTCTCTTCTTCCTTC";
        let motif = "TTCTT";
        let aligner = MotifAligner::align(target, motif);

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "{} {} {} {:.3} {}",
            bt.0,
            bt.1,
            bt.2,
            bt.3,
            underline_span(target, bt.1, bt.2)
        );

        assert_eq!(score, 132);
        assert_eq!(bt.1, 0);
        assert_eq!(bt.2, 84);
    }

    #[test]
    fn test_motif_alignment_real_two() {
        let target = "CCCTGTTTTTTTTTTTTTTTTTTTTTGGTAT";
        let motif = "T";

        let aligner = MotifAligner::align(target, motif);

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "motif:{} score:{} tstart:{} tend:{} pid:{:.3} tannoated:{}",
            motif,
            bt.0,
            bt.1,
            bt.2,
            bt.3,
            underline_span(target, bt.1, bt.2)
        );

        //aligner.print_matrix();
        assert_eq!(score, 42);
        assert_eq!(bt.1, 5);
        assert_eq!(bt.2, 26);
        assert_eq!(bt.3, 100.0);
    }

    #[test]
    fn test_motif_alignment_real_three() {
        let target = "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGAAAAACATTTATTCTGGGAATCTAACATCTGGACAAAAAATATTGGACTTTAGAAAAAGCCTACACCAAATTGTCTTATTTCTTCCCTCATACAGAATTCTGTACAGAATACATACAGAATAAGGAGGTGAAAAAAACTTTAAAAATAACATTGCTCCAGTTTGTCTGCAGGTATGTGATTTAAAATATCTCTGTTCTATTGAGGTATAGGCTGCAAACTTCGGTAAAATTAGAAAAATTAACAAACCCTTTTTAAAAAATACAAAACAAGTCCCCCCCTCAAAAAAGTCAAACAAACAAACAAAAAAACCTAACGTATCTTTCTTGGTACATATACACTAACTACAAAAAAATTTCATTTGGAAAATGTGTTCATTAAAAAAAAAAAAATCCCTCTTGTAAGTTCTCAATTGCCGAGTTAAAAAACGGCATTAAAATGATCTGTTAAACTCAAACTGTTCCCAAGCACCCGCCCCATGTCAGGAAGAACTGGATGACCGCCATTTCTAGAACCTGCTGCAGGGAAGCCAGCGGAGCAGGCGAGGTGGGGCTCCCAGAGGGTCGGCGGGCGCAGCTGCCTTGCGCGGGCTCCTAGGCGAGGGAGAGCCCCGCGCTGGGTAGGGGAAGTTCCTGCCCCCTCCGCTGCCCCACCCGCTAGGCCCCAGGCCCCAGATGTCTCCGTGCCGCCCTGTCTCTACCTGCGCTCTCTGCCCCGGTCCGCGAGGTCTAGGGTCACTACATCCCCTGTCTCGCGGGTACTCAGGGCTCCCTCTGCACTTAGAGCACAGTCACCACCCTGGACCGGTGTGGGTAGGCCCCCGCCTCCCCTCCTCGTCTGCTTCTACATGTCCACCCTTCTTGTGCCCGGCCTGTCCCCACCTTGTGTGCCTGTATGGCCGTGCTCAGTCCCGACCTCGGAAAGGCTGCTCCTTCTTGCGGTTGAGATCCTAGGACCAACACCCTCTCCTCACTCCATCCAGACTAGCGAATCCTCCCGCCCAAACCCCAGGTCTCACAGCCCCAGCACCAGGGAAAAGCATCTCATCAGGCCCTGCTGCGACAAGGAGAGCTCCGCGGGGCAGGGCCTTGCCTGTCTGGGCCCCAGCAAGGAACAAGTGCAAAGCTGAATCTGTAAAGTAATGGAACCCTAAAGGTGCTACCTGCTTCCCAGGGGAAGTCGTGAGAACTGATCTGTGTGCCTGTGTCCTACTGCCTTTAGTCAAGCGCCAGGCACACAGGAAAGAAACACTTGCTGCAAAAAAGCTACCGGCGCATCATTCATGCGGCCAGAGCCTTCACATAACAAATGCTTATTGCACACATCCTGTATTCCAGGCACCGTGTTATGAGTTTGGGATGGAGTACTCATGAAGACAGAATGATCCCTTTTGTGGAGCTTATAATGTCATGGGGAAAACATTAAATAAGTAAATACAAATGTTGATTAAAGTTGTCAGAAAACAAGTACACAGCATGCCACGGAAGCCTAGAGCAGCAGGACCCAACCTGTTGAAACATTTGTGTATGGGATCCTCAAATCCCAAATCCTGAGGCCTGCCTGTTCAGCCATGTAAACAGCCCCACCAGGAGGCCTGCCTTCCTTGCAGGTCTGTCAGGAGGCTCCTCTATTTTTGCCCAGCGCAGTCCCTCTCCAGGGCTCCCTCTTCAAAGCTAATGCTGGCTTCACTCCTCTCCTATACCCAAGAGATTTTAACTCAAAAGTGAAAACTTCTCCCAATCTGCCCCTTCTCCCTTGCTTCTAGGTACTTTCTAGATCTCCGAAGTTTTCTGATTTTTTCTAGACCACCACCCTCTACCTCCCCAGGGATTCCTTTTTATTTTTAAAAACTCTAAAGCTCCAGGCTCCTACAGCTGTGTGTGCACAAAAACAAACAAATAAAACCCTTAAGATTCTCCTTCTTCCCCACACATCTATTTCAAGGAAAGATGAGGGGGCAATTTTCCAATAGCACAAATCTCAGAATTCTGATAAATACAAATTAGAGATAGGGAAAAATACTATGATATACAAGCTATGATTTTAATTGGTATGAAACTTTGTGCAGAACTTTGACGATCATGGCTTGGAGAACTAGGAACTTTTCCTGTGTTCAGCAGGGACTTGGAAGAATTGCCCCCTTTCTATTACACAAAGATGCTTGCTAAAATATAGTGGTGGTTTGGGACCAAAGCAAGCTATGAGCACAATGACAATTATAGGAAAACAAGTGTGCCTCTGTACTGAAAACCCCAAGTCCAAGAATTTATTAGATTACCATCTTCTAAACATACTCTCCAGTAGCTCATGTATCTCTCACACCCGGTAAAAACTTAAGTTTCTGGGAAAAAAAGAATGTTTCTGAACATTAAATCATCATTTCAGGAACTTGCAACCTTAGTAAAGTCATGCTAAGAGTGAGTGTTTCCATGTG";
        let motif = "TG";
        let aligner = MotifAligner::align(target, motif);

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "motif:{} score:{} tstart:{} tend:{} tannoated:{}",
            motif,
            bt.0,
            bt.1,
            bt.2,
            underline_span(target, bt.1, bt.2)
        );

        aligner.print_matrix();
        assert_eq!(score, 56);
        assert_eq!(bt.1, 1);
        assert_eq!(bt.2, 29);
    }

    #[test]
    fn test_motif_alignment_real_four() {
        let target = "ATCGAATGGAATGGAATGGAATGGAATGGAATCAGCATAGAATGGAATCAAAAGCAATCATTCAATGGAATCTAATAGAATCATCGAATAGACTTGAATGTAATCATCATCGAATGGAGAAGAATGGAATCATCCAATGGACAGGAATGGAATCATCCTCGAATGGAATCGAATGGAATCATATAATGGACCCGAATGGAATCATCATTGAATGCAATAGAATGGAATCTTCATCGAATGGAATCATAAGGAATCATCTAATGGACACGAATGGAATCATCATCGAATTGAATAGAATGCAATCATCATCAAATGGAATCGAATGGAATCATCATCAAATAGAATTGAATGGAATCATCATGGAATGGAATTGAATGGAATCATCGAATGGACTCGAAAGCAATTATGCTCGAATGGAATCTAATGGAATCATCAAATGGACTCAAATGGAATCATCATCGAATGAAATCATATGGAATCATCGAATGAAACTGAATGGAATCATTGAATGGACTTGAAAGGAATTACTATGGAATGGGATTGAATGGAATCATTGAATGTACTTGAAAGGAATCATCATCAAGTGGAATCGAATGGAATCGTTGAATGGACTCGAATTGAATCTTTGAATGGAATCGAATGGAATCATCATTGAATGGAGTCAAATGGAATCATCATCAAATGGAATTGAATGGAATCTTCATTGAATGGACTCGAATGGAATCATCATCAAATGGAATCTAATCGAATCATCAATGAAGGGAATCGAATGGAATCATCATCGAATGGAATCGAATGGAATCATCAACAAGTGGAAACGAATGAAATCATCGAATGGAATCCAATGGTGTCATCGAATGGACACGAAAGGAATCATCGAATGGAATCACATGGAATCACCATCGAATGGAATCAAATGGAATCACCATCGAGTGGAATGTTATGAAATCATCTCATGGACTCGAAGGGAATCATCATCGAATGGAATCGAATGGAATCATTGAATGAAATGAAAGGAATCACCATCGAATGGAATGTTATGGAATGTTCTAATGGACTCGAAGGTAATCATCATCGAATGGAAACGAATGGAATCATTGAATGCAATTGAATGGAATCATCATCAAATGGAATCGAATGGAATCATCGAATGGAATCTGAATGGAATCATCAATGAATGGAATCAAATGGAATCATCTAATGGACGTGAATGGAATCATCATCGAAAGCAATGAAATGTAATTCAATCGAATGGACATGAATGGAATCATCATTGAATGGAATCAAATGGAATCCTCATCGAATGGAATCGAATGGAATCATCAAATGGAATAGAATGGAGTCATCGTCGAATGTAATCGAATGTAATCATCGAATGTCATCGAATGGAATCATCGACTGGAAAAGAATGGAATCATCATCGAATGGAAATGAATAGAATCACAGAATGAAATCGAATGGAATCATCATCGAATGGAGTCTAATGGAATAATAATCGAATGGAATAGAATGGAATCATCGAGTGGACACGAATGGAATCATCATTGAATGGATTCGAATAGAATCTTTTAATGAAATTGAATGGAATCAGCATGAAATGGAATCTAAAGGAATCATAGAATGGTATCGAATGGAATTATCATCGAATGGAATGGAATGGAATGGAATGGAATGGA";
        let motif = "GAATG";
        let aligner = MotifAligner::align(target, motif);

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "motif:{} score:{} tstart:{} tend:{} tannoated:{}",
            motif,
            bt.0,
            bt.1,
            bt.2,
            underline_span(target, bt.1, bt.2)
        );

        assert_eq!(score, 1425);
        assert_eq!(bt.1, 0);
        assert_eq!(bt.2, 1676);
        assert_float_relative_eq!(bt.3, 71.89, 0.01);
    }

    #[test]
    fn test_motif_alignment_real_five() {
        let target = "AGAAAGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGATGGGAGTGTCTGTGTATGTGTGTGAGAGAGAGAGATGGGAGACTGTGTGTGTGTACGCGAGAGAGATGGGAGAGTGTGTGTGAGAGAGATGGGGAGAGTGTGTGTGTGTGTACGAGAGAGATG";
        let motif = "GT";
        let aligner = MotifAligner::align(target, motif);

        let (score, _i, _j) = aligner.best_score();

        let bt = aligner.backtrace();

        println!(
            "motif:{} score:{} tstart:{} tend:{} tannoated:{}",
            motif,
            bt.0,
            bt.1,
            bt.2,
            underline_span(target, bt.1, bt.2)
        );

        aligner.print_matrix();
        assert_eq!(score, 143);
        assert_eq!(bt.1, 5);
        assert_eq!(bt.2, 148);
    }

    #[test]
    fn test_mismatch_penalty() {
        let target = "AAAA";
        let motif = "TTTT";
        let aligner = MotifAligner::align(target, motif);

        let (score, _, _) = aligner.best_score();

        // No matches, all mismatches: best local alignment should be 0 (Smith-Waterman)
        assert_eq!(score, 0);
    }

    #[test]
    fn test_perfect_match() {
        let target = "ACGT";
        let motif = "ACGT";
        let aligner = MotifAligner::align(target, motif);
        let (score, _, _) = aligner.best_score();
        assert_eq!(score, 8); // 4 matches × 2
    }

    #[test]
    fn test_perfect_match_2() {
        let target = "CCCATATATATATATATATATATATCCCC";
        let motif = "AT";
        let aligner = MotifAligner::align(target, motif);
        let (score, _, _) = aligner.best_score();
        let bt = aligner.backtrace();
        aligner.print_matrix();
        println!("{}", underline_span(target, bt.1, bt.2));
        assert_eq!(score, 44);
        assert_eq!(bt.3, 100.0);
        assert_eq!(bt.1, 3);
        assert_eq!(bt.2, 25);
    }
}
