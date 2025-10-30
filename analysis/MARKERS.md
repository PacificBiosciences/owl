# Guide to build marker set

1. Pull UCSC genome browser repeat masker annotations

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
gunzip rmsk.txt.gz
```

2. Select markers using python script

```
python3 ./scripts/repfilter.py \
  --min-length 20 \
  --max-length 100 \
  --min-identity 99.99 \
  --max-unit-len 5 \
  rmsk.txt \
  > repeats.filtered.bed
``` 

3. Get main chromosomes

```
grep -v -E "Un|random|alt|fix|chrY" repeats.filtered.bed > repeats.mainchroms.bed
```


4. Subtract “odd” regions (e.g., assembly gaps or masked regions)

```
bedtools subtract -A \
  -a repeats.mainchroms.bed \
  -b ./data/GRCh38.oddRegions.bed.gz \
  > repeats.noOdd.bed
```

5. S