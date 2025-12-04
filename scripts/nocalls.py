#!/usr/bin/env python3
import sys
import argparse


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Emit rows where ALL haplotype entries have missing CV ('.') "
            "for Owl v0.3-style tables."
        )
    )
    p.add_argument(
        "infile",
        nargs="?",
        default="-",
        help="Input TSV (default: stdin)",
    )
    p.add_argument(
        "--print-header",
        action="store_true",
        help="Print header lines (starting with '#')",
    )
    p.add_argument(
        "--summary",
        action="store_true",
        help="Print summary at the end to stderr",
    )
    return p.parse_args()


def open_in(path):
    if path == "-" or path is None:
        return sys.stdin
    return open(path, "r")


def all_cv_missing_in_sample(sample_str, cv_pos):
    """
    Given one sample column from the Owl file and the index of 'CV' in the
    format string, return True if ALL haplotypes have CV == '.'.

    Sample format:
      PS,HP,CT,MU,CV,LV;PS,HP,CT,MU,CV,LV;...

    Examples:
      '0,0,0,.,.,.'                  -> all missing
      '257717,0,6,44.5,21.5,.'       -> not missing (CV=21.5)
    """
    sample_str = sample_str.strip()
    if sample_str == "" or sample_str == ".":
        # No data at all; treat as missing
        return True

    # Semicolon-separated haplotype entries
    for entry in sample_str.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(",")
        if cv_pos >= len(parts):
            # No CV field in this hap entry – treat as missing for this hap
            continue
        cv_val = parts[cv_pos].strip()
        if cv_val != "." and cv_val != "":
            # At least one haplotype has a real CV value
            return False

    # All haplotypes had '.' (or no CV field)
    return True


def process_stream(fh, print_header=False, summary=False):
    kept = 0
    total = 0

    for raw in fh:
        line = raw.rstrip("\n")

        # Header / meta lines
        if line.startswith("#"):
            if print_header:
                print(line)
            continue

        if not line.strip():
            continue

        fields = line.split("\t")
        if len(fields) < 4:
            # Need at least Region, info, format, and one sample
            continue

        format_str = fields[2]
        fmt_keys = format_str.split(":")

        try:
            cv_pos = fmt_keys.index("CV")
        except ValueError:
            # No CV field in this row; skip from consideration
            continue

        total += 1

        # Sample columns start at index 3
        all_missing = True
        for sample_str in fields[3:]:
            if not all_cv_missing_in_sample(sample_str, cv_pos):
                all_missing = False
                break

        if all_missing:
            print(line)
            kept += 1

    if summary:
        print(
            f"# kept {kept} / {total} rows where ALL haplotype CV values are missing ('.')",
            file=sys.stderr,
        )


def main():
    args = parse_args()
    with open_in(args.infile) as fh:
        process_stream(
            fh,
            print_header=args.print_header,
            summary=args.summary,
        )


if __name__ == "__main__":
    main()
