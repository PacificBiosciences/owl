#!/usr/bin/env python3
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Emit rows where ALL haplotype entries have missing CV ('.')."
    )
    p.add_argument("infile", nargs="?", default="-", help="Input TSV (default: stdin)")
    p.add_argument("--print-header", action="store_true",
                   help="Print header line first if present")
    p.add_argument("--summary", action="store_true",
                   help="Print summary at the end to stderr")
    return p.parse_args()

def open_in(path):
    return sys.stdin if path == "-" else open(path, "r", encoding="utf-8")

def all_cv_missing(sample_str, cv_idx):
    """
    Return True iff every haplotype entry has CV missing ('.' or empty).
    - Missing/empty sample column counts as missing for all.
    - If an entry is too short to contain CV, treat it as missing (fail).
    """
    if not sample_str or sample_str == ".":
        return True
    for entry in sample_str.split(";"):
        entry = entry.strip()
        if not entry:
            # empty entry -> missing CV
            continue
        fields = [f.strip() for f in entry.split(",")]
        if cv_idx >= len(fields):
            # no CV position => treat as missing (fail)
            continue
        v = fields[cv_idx]
        # If CV present and not '.', then this entry passes -> row shouldn't be emitted
        if v not in ("", "."):
            return False
    return True

def main():
    args = parse_args()
    total = kept = 0
    header_emitted = False
    region_i = fmt_i = sample_i = None

    with open_in(args.infile) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("##"):
                continue

            if line.startswith("#"):
                cols = line[1:].split("\t")
                # column indices (best-effort fallbacks)
                region_i = cols.index("Region") if "Region" in cols else 0
                fmt_i    = cols.index("format") if "format" in cols else 3
                sample_i = len(cols) - 1  # single sample column is last
                if args.print_header:
                    print(line[1:])
                    header_emitted = True
                continue

            cols = line.split("\t")
            if fmt_i is None or sample_i is None:
                # no header case: assume Region, len, motif, format, sample
                if len(cols) < 5:
                    continue
                region_i, fmt_i, sample_i = 0, 3, 4

            fmt = cols[fmt_i]
            sample_str = cols[sample_i]

            # find CV index from format, e.g. "hp:ct:mu:cv:ln"
            fmt_parts = [f.strip().lower() for f in fmt.split(":")]
            try:
                cv_idx = fmt_parts.index("cv")
            except ValueError:
                # no CV in format: treat as not evaluable -> skip row
                continue

            total += 1
            if all_cv_missing(sample_str, cv_idx):
                print(line)
                kept += 1

    if args.summary:
        print(f"# kept {kept} / {total} rows where ALL haplotype CV values are missing ('.')",
              file=sys.stderr)

if __name__ == "__main__":
    main()
