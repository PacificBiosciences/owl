#!/usr/bin/env python3

import csv
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Process motif table.")
    parser.add_argument(
        "-i", "--input",
        metavar="FILE",
        type=str,
        required=True,
        help="Input tab-delimited file from owl"
    )
    args = parser.parse_args()

    # Read all lines so we can reliably take the *last* header line
    with open(args.input, newline='') as fh:
        lines = fh.readlines()

    # Find the last header line beginning with '#'
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('#'):
            header_idx = i
    if header_idx is None:
        print("Error: no header line starting with '#' found.", file=sys.stderr)
        sys.exit(1)

    # Parse header fields (strip leading '#' and split on tabs if present, else whitespace)
    raw_header = lines[header_idx].lstrip('#').strip()
    if '\t' in raw_header:
        header = raw_header.split('\t')
    else:
        header = raw_header.split()

    # Validate required columns
    try:
        length_col = header.index("len")
        motif_col = header.index("motif")
        format_col = header.index("format")
    except ValueError as e:
        print(f"Error: required column missing in header: {e}", file=sys.stderr)
        sys.exit(1)

    # Keep ONLY the last header entry (last column), which should be the sample name
    if len(header) <= format_col + 1:
        print("Error: no sample columns found after 'format'.", file=sys.stderr)
        sys.exit(1)
    sample_idx = len(header) - 1
    sample_name = header[sample_idx]

    # Create a CSV reader over the data lines following the header
    data_lines = lines[header_idx + 1 :]
    reader = csv.reader(data_lines, delimiter="\t")

    for row in reader:
        if not row or row[0].startswith('#'):
            continue  # skip any stray headers/blank lines

        # Defensive: ensure row has enough columns
        if len(row) <= sample_idx:
            continue

        region = row[0]
        length = row[length_col]
        motif  = row[motif_col]
        sample_data = row[sample_idx]

        # Split by ';' to get each haplotype entry
        entries = [e for e in sample_data.split(";") if e.strip()]
        for e in entries:
            fields = e.split(",")
            # Expect at least: hp, ct, mu, cv
            if len(fields) >= 4:
                hp, ct, mu, cv = fields[0:4]

                # Skip if any is '.' or hp == '0'
                if "." in (length, motif, hp, mu, cv):
                    continue

                print(f"{region}\t{length}\t{motif}\t{sample_name}\t{ct}\t{hp}\t{mu}\t{cv}")

if __name__ == "__main__":
    main()
