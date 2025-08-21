# Owl

<h1 align="center"><img width="300px" src="logo/owl-logo.svg"/></h1>

<h1 align="center">Microsatellite instability (MSI) analysis for HiFi data.</h1>



Authors: [Zev Kronenberg](https://github.com/zeeev), [Khi Pin Chua](https://github.com/proteinosome), [Egor Dolzhenko](https://github.com/egor-dolzhenko), [Mark Chaisson](https://github.com/mchaisso) 

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
After running `owl score` there are two output files, `{prefix}.owl-motif-counts.txt` and `{prefix}.owl-scores.txt`. The scores file provides a summary score for each sample, whereas the motif file breaks down the score by motif. 

### Summary score output
| sample  | #high | #low   | %high |
|---------|------|-------|-------|
| NA12878 | 1868 | 72231 | 2.521 |
|   ...   |  ... | ...   | ...   |

The high and low counts are the number of loci that have a high coefficient of variance, the main measure of MSI.

### Summary of motif output
| motif  | #high | #low | %high |
|--------|------|-----|-------|
| TAGGAC | 0    | 0   | 0.00  |
| ... | ...    | ...   | ...  |

The motif file contains the same information, but summarizes the same information across motifs. If multiple samples are scored together, the motif stats 


## Changelog 

* v0.1.1 -- August 21
  - initial release to github

## Questions, Comments, Feedback:
Please feel free to file a ticket at:
[PacBio GitHub Issues](https://github.com/PacificBiosciences/pbbioconda/issues)