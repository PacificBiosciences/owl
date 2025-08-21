# Owl

<h1 align="center"><img width="300px" src="logo/owl-logo.svg"/></h1>

<h1 align="center">Microsatellite instability analysis for HiFi data.</h1>



Authors: [Zev Kronenberg](https://github.com/zeeev), [Khi Pin Chua](https://github.com/proteinosome), [Egor Dolzhenko](https://github.com/egor-dolzhenko), [Mark Chaisson](https://github.com/mchaisso) 


## Running

### step one, profile the repeats
```
gunzip data/Simple-repeats-50k.filt.bed.gz

 owl profile --bam NA12878.haplotagged.bam --regions data/Simple-repeats-50k.bed --sample NA12878 > NA12878.results.txt
```

### step two, summarize and score sample
```
owl score --file NA12878.results.txt --prefix NA12878-scored

```

## Output of step one (profiling)

| Column            | Details                                                                 |
|-------------------|-------------------------------------------------------------------------|
| Region            | Simple region format: chrom:start-end                                   |
| Region length     | Repeat masker annotation length in reference.                           |
| Motif             | Motif to search for in the reads                                        |
| Format            | Format information per haplotype label                                  |

## Format information
| Field             | key | Details                                       |
|-------------------|-----|-----------------------------------------------|
| Haplotype         | hp  | Haplotype label from bam file (0,1,2,...)     |
| Read Count        | ct  | Number of reads after filtering               |
| Mean Length       | mu  | Mean length of motif span across reads        |
| Coefficient of variation | cv  | Coefficient of variation of motif span              | 
| Repeat lengths    | ln  | Colon separated list of motif spans           |
       

## Building target regions

```
 #Download repeat masker file for GRCh38
 wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
 #Decompress 
 gunzip rmsk.txt.gz

 #Keeping simple repeats with a divergence less than 50 (5%), > 20 bp && < 300 bp, and motifs smaller than 10.
 grep Simple rmsk.txt |  perl -lane 'print if $F[2] <  50 && $F[7] - $F[6] > 20 ' | perl -lane 'print if $F[7] - $F[6] < 300' | perl -lane  'print "$F[5]\t$F[6]\t$F[7]\t$F[10]"' | perl -lane 'print if length($F[-1]) <= 9' | sort -k1,1 -k2,2n  >  Simple-repeats.bed
 perl -lane 'print "$F[5]\t$F[6]\t$F[7]\t$F[10]"' rmsk.txt | sort -k1,1 -k2,2n > all-reps.bed

 #Merge overlapping annotations
 bedtools merge -i all-reps.bed -c 4 -o distinct | grep "," | bedtools sort -i - > all-ovl-reps.bed

 #Removing repeats with close neighbors, filtering chromosomes, and sampling 50k sites.
 bedtools closest -d -io -a Simple-repeats.bed -b all-reps.bed | perl -lane 'print if $F[-1] > 100' | bedtools subtract -A -a - -b all-ovl-reps.bed   | grep -v 'fix' | grep -v alt |  grep -v chrUn | grep -v random | grep -v chrY  | shuf | head -n 50000 | bedtools sort -i - | cut -f 1-4 | uniq > Simple-repeats-50k.bed

 Removing segmental duplications
 bedtools subtract -A -a Simple-repeats-50k.bed -b GRCh38.segdups.bed.gz > Simple-repeats-50k.filt.bed

```
