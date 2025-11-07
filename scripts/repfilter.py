#!/usr/bin/env python3
"""
rmsk_parse.py — simple parser & filter for UCSC RepeatMasker (rmsk.txt[.gz])

Input format (17 tab-separated columns):
  1 bin, 2 swScore, 3 milliDiv, 4 milliDel, 5 milliIns,
  6 genoName, 7 genoStart, 8 genoEnd, 9 genoLeft,
 10 strand, 11 repName, 12 repClass, 13 repFamily,
 14 repStart, 15 repEnd, 16 repLeft, 17 id
"""
from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from dataclasses import dataclass, asdict
from typing import Iterable, Iterator, Optional, TextIO, List

COLS = [
    "bin","swScore","milliDiv","milliDel","milliIns",
    "genoName","genoStart","genoEnd","genoLeft",
    "strand","repName","repClass","repFamily",
    "repStart","repEnd","repLeft","id"
]

@dataclass
class RMSKRecord:
    bin: int
    swScore: int
    milliDiv: int             # divergence in thousandths of a percent
    milliDel: int
    milliIns: int
    genoName: str
    genoStart: int
    genoEnd: int
    genoLeft: int
    strand: str
    repName: str
    repClass: str
    repFamily: str
    repStart: int
    repEnd: int
    repLeft: int
    id: int

    @property
    def length(self) -> int:
        return self.genoEnd - self.genoStart

    @property
    def motif(self) -> Optional[str]:
        m = re.fullmatch(r"\(([ACGTN]+)\)n", self.repName)
        return m.group(0) if m else None

    @property
    def motif_core(self) -> Optional[str]:
        m = re.fullmatch(r"\(([ACGTN]+)\)n", self.repName)
        return m.group(1) if m else None

    @property
    def identity(self) -> float:
        # milliDiv is "divergence * 1000"
        # identity (%) = 100 - divergence(%)
        #              = 100 - (milliDiv / 10)
        return 100.0 - (self.milliDiv / 10.0)

def open_maybe_gzip(path: str) -> TextIO:
    if path == "-" or path is None:
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "rt")

def parse_rmsk(handle: TextIO) -> Iterator[RMSKRecord]:
    for line in handle:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split()
        if len(parts) != 17:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 17:
                continue
        try:
            yield RMSKRecord(
                bin=int(parts[0]),
                swScore=int(parts[1]),
                milliDiv=int(parts[2]),
                milliDel=int(parts[3]),
                milliIns=int(parts[4]),
                genoName=parts[5],
                genoStart=int(parts[6]),
                genoEnd=int(parts[7]),
                genoLeft=int(parts[8]),
                strand=parts[9],
                repName=parts[10],
                repClass=parts[11],
                repFamily=parts[12],
                repStart=int(parts[13]),
                repEnd=int(parts[14]),
                repLeft=int(parts[15]),
                id=int(parts[16]),
            )
        except ValueError:
            continue

def is_basic_motif(core: str) -> bool:
    """
    Keep if:
      - motif length ≤ 6
      - AND one of:
          * homopolymer (len==1, A/C/G/T)
          * dinucleotide (len==2, A/C/G/T)
          * contains any run of the same base like 'AA','CC','GG','TT'
            (e.g., ATT, GCC, AATA)
    """
    if not core:
        return False
    core = core.upper()
    if not re.fullmatch(r"[ACGT]+", core):
        return False
    if len(core) > 6:
        return False
    if len(core) == 1:
        return True
    if len(core) == 2:
        return True
    return re.search(r"(.)\1", core) is not None

def filter_records(
    records: Iterable[RMSKRecord],
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
    rep_class: Optional[str],
    rep_family: Optional[str],
    motif_regex: Optional[str],
    min_len: Optional[int],
    max_len: Optional[int],
    min_identity: Optional[float],
    motif_core_basic: bool,
    max_unit_len: Optional[int],
) -> Iterator[RMSKRecord]:
    motif_pat = re.compile(motif_regex) if motif_regex else None
    for r in records:
        if chrom and r.genoName != chrom:
            continue
        if start is not None and r.genoEnd <= start:
            continue
        if end is not None and r.genoStart >= end:
            continue
        if rep_class and r.repClass != rep_class:
            continue
        if rep_family and r.repFamily != rep_family:
            continue
        if min_len is not None and r.length < min_len:
            continue
        if max_len is not None and r.length > max_len:
            continue

        # Motif-based filters
        if motif_pat or motif_core_basic or (max_unit_len is not None):
            core = r.motif_core
            if core is None:
                continue  # require "(...)n"
            if motif_pat and not motif_pat.search(r.motif or r.repName):
                continue
            if motif_core_basic and not is_basic_motif(core):
                continue
            if (max_unit_len is not None) and (len(core) > max_unit_len):
                continue

        # Identity filter
        if min_identity is not None:
            ident = r.identity  # use the property (100 - milliDiv/1000)
            if ident < min_identity:
                continue

        yield r

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Parse and filter UCSC RepeatMasker rmsk.txt(.gz) and print records."
    )
    p.add_argument("input", help="Path to rmsk.txt or rmsk.txt.gz (use '-' for stdin)")
    p.add_argument("--chrom", help="Filter: chromosome (e.g. chr1)")
    p.add_argument("--start", type=int, help="Filter: start (0-based, inclusive) on genome")
    p.add_argument("--end", type=int, help="Filter: end (0-based, exclusive) on genome")
    p.add_argument("--class", dest="rep_class", help="Filter: repClass (e.g. Simple_repeat, SINE, LTR)")
    p.add_argument("--family", dest="rep_family", help="Filter: repFamily (e.g. Alu, ERVL-MaLR)")
    p.add_argument("--motif-regex", help="Filter regex applied to motif like '(AC)n' if present, else repName")
    p.add_argument(
        "--motif-core-basic",
        action="store_true",
        help=(
            "Keep only Simple_repeat motifs that (1) have length ≤ 6 and "
            "(2) are homopolymers, dinucleotides, or contain any run of the "
            "same base (e.g., ATT, GCC, AATA)."
        ),
    )
    p.add_argument(
        "--max-unit-len",
        type=int,
        default=None,
        help=(
            "Keep only Simple_repeat motifs of the form '(XXX)n' whose repeat "
            "unit length is ≤ this value (e.g., 1 for homopolymers, 2 for di-, …)."
        ),
    )
    p.add_argument("--min-length", type=int, default=None, help="Filter: minimum genomic length")
    p.add_argument("--max-length", type=int, default=None, help="Filter: maximum genomic length")
    p.add_argument(
        "--min-identity",
        type=float,
        help="Minimum percent identity (0–100). Computed as 100 - milliDiv/1000.",
    )
    p.add_argument(
        "--fields",
        default="genoName,genoStart,genoEnd,repName,repClass,repFamily",
        help=(
            "Comma-separated fields to print (default: %(default)s). "
            "Available: "
            + ",".join(COLS)
            + " + length,motif,motif_core,identity. "
            "Regardless of this setting, an 'identity' column (percent) "
            "is appended as the final column."
        ),
    )
    args = p.parse_args(argv)

    with open_maybe_gzip(args.input) as fh:
        recs = parse_rmsk(fh)
        recs = filter_records(
            recs,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
            rep_class=args.rep_class,
            rep_family=args.rep_family,
            motif_regex=args.motif_regex,
            min_len=args.min_length,
            max_len=args.max_length,
            min_identity=args.min_identity,
            motif_core_basic=args.motif_core_basic,
            max_unit_len=args.max_unit_len,
        )

        wanted = [f.strip() for f in args.fields.split(",") if f.strip()]
        for r in recs:
            row = []
            d = asdict(r)

            # user-requested fields
            for f in wanted:
                if f in ("length", "motif", "motif_core", "identity"):
                    val = getattr(r, f)
                elif f in d:
                    val = d[f]
                else:
                    sys.stderr.write(f"[warn] unknown field '{f}'\n")
                    val = ""
                row.append(f"{val:.3f}" if isinstance(val, float) else str(val if val is not None else ""))

            # always append identity as last column
            row.append(f"{r.identity:.3f}")

            print("\t".join(row))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
