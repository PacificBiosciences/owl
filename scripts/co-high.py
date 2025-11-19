#!/usr/bin/env python3

import sys
from collections import defaultdict

# Fixed column indices based on your file layout:
# chr  refsize  motif  sample  coverage  haplotype  mean  cv  number-haplotypes
COL_SITE  = 0   # marker
COL_MOTIF = 2   # motif
COL_HAP   = 5   # haplotype
COL_CV    = 7   # cv

THRESHOLD = 5.0
HAPS = ["0", "1", "2"]


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} <input_file> <sample_name>\n")
        sys.exit(1)


    sample = sys.arv[2]
    infile = sys.argv[1]

    # site -> set of high haplotypes
    high_by_site = defaultdict(set)
    # site -> motif
    motif_by_site = {}

    # ---- Parse input file and collect high haplotypes ----
    with open(infile) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            # Expect at least 9 columns
            if len(parts) < 9:
                continue

            site   = parts[COL_SITE]
            motif  = parts[COL_MOTIF]
            hap    = parts[COL_HAP]
            cv_str = parts[COL_CV]

            motif_by_site[site] = motif

            if cv_str == ".":
                continue

            try:
                cv_val = float(cv_str)
            except ValueError:
                # Header or malformed line; skip
                continue

            if cv_val > THRESHOLD:
                high_by_site[site].add(hap)

    # ---- 3x3 matrix for haplotypes 0,1,2 ----
    matrix = {h1: {h2: 0 for h2 in HAPS} for h1 in HAPS}

    # Track co-high sites: (site, motif, [high_haps])
    cohigh_sites = []

    # Consider all sites we know about (union of those with motif or highs)
    all_sites = sorted(set(high_by_site.keys()) | set(motif_by_site.keys()))

    for site in all_sites:
        highs = high_by_site.get(site, set()).intersection(HAPS)
        if not highs:
            continue

        highs_list = sorted(highs)

        # Diagonal: each high hap gets +1
        for h in highs_list:
            matrix[h][h] += 1

        # Off-diagonal: if >= 2 highs, track site and increment pair counts
        if len(highs_list) >= 2:
            motif = motif_by_site.get(site, ".")
            cohigh_sites.append((site, motif, highs_list))

            for i in range(len(highs_list)):
                for j in range(i + 1, len(highs_list)):
                    hi = highs_list[i]
                    hj = highs_list[j]
                    matrix[hi][hj] += 1
                    matrix[hj][hi] += 1

    # ---- Print 3x3 matrix ----
    print("hap\t" + "\t".join(HAPS))
    for h1 in HAPS:
        row = [h1] + [str(matrix[h1][h2]) for h2 in HAPS]
        print("\t".join(row))

    # ---- Print co-high sites with motif ----
    print("\n# sites_with_cohigh")
    for site, motif, haps in cohigh_sites:
        print(f"{site}\t{motif}\t{','.join(haps)}")

    # ---- Fraction for haplotypes 1 and 2 ----
    h1 = "1"
    h2 = "2"

    sites_h1 = 0
    sites_h2 = 0
    sites_h1_or_h2 = 0
    sites_h1_and_h2 = 0

    # Only iterate over sites where at least one hap was high (0,1,2)
    for site, highs in high_by_site.items():
        has1 = (h1 in highs)
        has2 = (h2 in highs)

        if has1:
            sites_h1 += 1
        if has2:
            sites_h2 += 1

        if has1 or has2:
            sites_h1_or_h2 += 1
        if has1 and has2:
            sites_h1_and_h2 += 1

    if sites_h1_or_h2 > 0:
        frac_both_given_either = sites_h1_and_h2 / sites_h1_or_h2
    else:
        frac_both_given_either = 0.0

    # Expected fraction under independence:
    # p1 = P(h1 high), p2 = P(h2 high) over all sites with any hap high
    n_sites_for_p = len(high_by_site)
    if n_sites_for_p > 0:
        p1 = sites_h1 / n_sites_for_p
        p2 = sites_h2 / n_sites_for_p
        denom = (p1 + p2 - p1 * p2)
        if denom > 0:
            expected_frac = (p1 * p2) / denom
        else:
            expected_frac = 0.0
    else:
        p1 = p2 = expected_frac = 0.0

    ratio = (frac_both_given_either / expected_frac) if expected_frac > 0 else 0.0

    print("\n# Stats: ")
    print(f"stats\t{sample}\tsites_high\t{n_sites_for_p}")
    print(f"stats\t{sample}\tsites_h1_high\t{sites_h1}")
    print(f"stats\t{sample}\tsites_h2_high\t{sites_h2}")
    print(f"stats\t{sample}\tunion\t{sites_h1_or_h2}")
    print(f"stats\t{sample}\tintersection\t{sites_h1_and_h2}")
    print(f"stats\t{sample}\tjaccard_index\t{frac_both_given_either:.6f}")
    print(f"stats\t{sample}\tp1 = P(h1 high)\t{p1:.6f}")
    print(f"stats\t{sample}\tp2 = P(h2 high)\t{p2:.6f}")
    print(f"stats\t{sample}\texpected_frac_independent\t{expected_frac:.6f}")
    print(f"stats\t{sample}\tobserved_over_expected\t{ratio:.6f}")


if __name__ == "__main__":
    main()
