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
    parser.add_argument(
        "--max",
        action="store_true",
        help=(
            "If set, for each row only output the haplotype entry with the "
            "maximum numeric CV."
        ),
    )
    args = parser.parse_args()

    # Read all lines so we can reliably take the *last* header line
    with open(args.input, newline="") as fh:
        lines = fh.readlines()

    # Find the last header line beginning with '#'
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("#"):
            header_idx = i
    if header_idx is None:
        print("Error: no header line starting with '#' found.", file=sys.stderr)
        sys.exit(1)

    # Parse header fields (strip leading '#' and split on tabs if present, else whitespace)
    raw_header = lines[header_idx].lstrip("#").strip()
    if "\t" in raw_header:
        header = raw_header.split("\t")
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

    # Assume the last column after 'format' is the sample data column
    if format_col >= len(header) - 1:
        print("Error: no sample columns found after 'format'.", file=sys.stderr)
        sys.exit(1)
    sample_idx = len(header) - 1
    sample_name = header[sample_idx]

    # Create a CSV reader over the data lines following the header
    data_lines = lines[header_idx + 1 :]
    reader = csv.reader(data_lines, delimiter="\t")

    for row in reader:
        if not row or row[0].startswith("#"):
            continue  # skip any stray headers/blank lines

        # Defensive: ensure row has enough columns
        if len(row) <= sample_idx:
            continue

        region = row[0]
        length = row[length_col]
        motif = row[motif_col]
        sample_data = row[sample_idx]

        # Split by ';' to get each haplotype entry
        entries = [e for e in sample_data.split(";") if e.strip()]
        if not entries:
            continue

        passing = 0

        # For --max mode
        best_fields = None   # (hp, ct, mu, cv_str)
        best_cv_val = None
        has_numeric_cv = False
        fallback_fields = None

        # First pass: compute passing count and, if needed, track max CV
        for e in entries:
            fields = e.split(",")
            # Expect at least: hp, ct, mu, cv
            if len(fields) < 4:
                continue
            hp, ct, mu, cv = fields[0:4]

            # First valid entry as a fallback if we never find a numeric CV
            if fallback_fields is None:
                fallback_fields = (hp, ct, mu, cv)

            if cv != ".":
                passing += 1
                try:
                    cv_val = float(cv)
                except ValueError:
                    cv_val = None

                if cv_val is not None:
                    if (not has_numeric_cv) or (cv_val > best_cv_val):
                        has_numeric_cv = True
                        best_cv_val = cv_val
                        best_fields = (hp, ct, mu, cv)

        if args.max:
            # Only one output per row
            if has_numeric_cv and best_fields is not None:
                hp, ct, mu, cv = best_fields
            elif fallback_fields is not None:
                hp, ct, mu, cv = fallback_fields
            else:
                continue

            print(
                f"{region}\t{length}\t{motif}\t{sample_name}"
                f"\t{ct}\t{hp}\t{mu}\t{cv}\t{passing}"
            )
        else:
            # Original behaviour: print one line per haplotype entry
            for e in entries:
                fields = e.split(",")
                # Expect at least: hp, ct, mu, cv
                if len(fields) < 4:
                    continue
                hp, ct, mu, cv = fields[0:4]
                print(
                    f"{region}\t{length}\t{motif}\t{sample_name}"
                    f"\t{ct}\t{hp}\t{mu}\t{cv}\t{passing}"
                )


if __name__ == "__main__":
    main()
