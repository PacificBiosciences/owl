<h1 align="center">Owl</h1>

<h1 align="center"><img width="300px" src="logo/owl-logo.svg"/></h1>

<h1 align="center">Microsatellite instability (MSI) analysis for HiFi data.</h1>



Authors: [Zev Kronenberg](https://github.com/zeeev), [Khi Pin Chua](https://github.com/proteinosome), [Egor Dolzhenko](https://github.com/egor-dolzhenko), [Mark Chaisson](https://github.com/mchaisso), [Mike Eberle]() 

## Building
```
 cd owl
 cargo build --release
```

## Running

### step one, profile the repeats
```
 # decompress default simple repeat catalog.
 gunzip data/Simple-repeats-50k.filt.bed.gz

 # profile the repeats.
 owl profile --bam NA12878.haplotagged.bam --regions data/Simple-repeats-50k.bed --sample NA12878 > NA12878.results.txt
```

### step two, summarize and score sample(s)
```
# run the scoring step.
owl score --file NA12878.results.txt --prefix NA12878-scored
```

## Output
After running `owl score` there are two output files, `{prefix}.owl-motif-counts.txt` and `{prefix}.owl-scores.txt`. The score file provides a summary MSI score for each sample, whereas the motif file breaks down the score by motif. 

| sample   | #high | #low | %high | #phased | %phased | #sites | #passing | %passing | qc   |
| -------- | ----: | ---: | ----: | ------: | ------: | -----: | -------: | -------: | :--- |
| NA12878  |    15 |  758 |  1.94 |     378 |   75.60 |    500 |      399 |    79.80 | pass |
| …        |     … |    … |     … |      …  |      …  |      … |        … |        … | …    |

High and low are counts of haplotypes (multiple per loci) with high vs. low coefficient of variation (CV). %high is the proportion of loci with high CV (our primary MSI metric). QC reflects data completeness: it reports the percentage of sites with reliable measurements (%passing), and the qc column labels each sample pass or fail based on that percentage.

### Summary of motif output
| motif  | #high | #low | %high |
|--------|------|-----|-------|
| TAGGAC | 0    | 0   | 0.00  |
| ... | ...    | ...   | ...  |

The motif file contains the same information, but summarizes the same information across motifs. If multiple samples are scored together, the motif stats 


## Changelog 
* v0.3.0 -- Sept 22 2025
  - Warn if `owl profile` does encounter a phased region.
  - Switch polarity of un-phased read filter, keep unphased reads unless region contains phased reads.
  - Update score report to include phasing information.
  - Fix panic on bam open.
* v0.2.1 -- Sept 3 2025
  - Fix duplicate header field
* v0.2.0 -- August 28 2025
  - Add QC metric to reporting
* v0.1.3 -- August 26 2025
  - fix memory reporting bug
* v0.1.2 -- August 21 2025
  - initial release to github

## Questions, Comments, Feedback:
Please feel free to file a ticket at:
[PacBio GitHub Issues](https://github.com/PacificBiosciences/pbbioconda/issues)