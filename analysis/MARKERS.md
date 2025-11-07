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
  --min-identity 97 \
  --max-unit-len 6 \
  rmsk.txt \
  > repeats.filtered.bed
``` 

3. Get main chromosomes

```
grep -v -E "Un|random|alt|fix|chrY" repeats.filtered.bed | cut -f 1-4 > repeats.mainchroms.bed
```

Following these steps, additional refinements can be applied to the marker set. For instance, one can evaluate the frequency of site no-calls or identify and remove markers exhibiting excessive variability across populations.