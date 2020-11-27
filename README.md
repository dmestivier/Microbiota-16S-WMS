This is our repository of the scripts and data 
we generated and create during the work
published in our manuscript:

```

The limits and avoidance of biases in metagenomic analyses of human fecal
microbiota 

Emma Bergsten, Denis Mestivier and Iradj Sobhani
```

## Raw data

### 16S

Directory: `./Tables/16S`

- QIIME1: otu_taxa_EC2M3_16S_16S_clean_silva.csv
- DADA2: otu_taxa_EC2M3_16S_DADA2-Filt150_Silva128_counts-and-annotation.csv
- MALT+MEGAN6: Comparison-EC2M3-16S_MALT-SILVA128-fromSAM_rev-19-01-2018_Taxonomy-mod.csv

These table contain abundances and taxonomy
(row: taxon / col: samples)

### 16S / Merge

REP: `./Tables/16S`

Data are merged and summarized at genus level in one table with or without
normalisation.

- `agregation-MALT-Qiime-DADA2_CCR1-16S.R`
- `agregation-norm-MALT-Qiime-DADA2_CCR1-16S.R`

```
CMP_CCR116S-Normaux-MALT-QIIME-DADA2_counts_rev-22-11-2019.csv
CMP_CCR116S-Normaux-MALT-QIIME-DADA2_taxa_rev-22-11-2019.csv
CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_counts_rev-22-11-2019.csv
CMP_CCR116S-Normaux-norm-MALT-QIIME-DADA2_taxa_rev-22-11-2019.csv
```

The file contains 89 Genus shared by the three pipelines.

---------------------------------------------------------------

### WGS

REP : `./Tables/WGS`

Tables x :

- DIAMOND+MEGAN6: Comparison-EC2M3-WGS_Taxonomy-R1_rev-17-05-2017_bc-taxo.csv
- DIAMOND+MEGAN6: Comparison-EC2M3-WGS_Taxonomy-R1_rev-17-05-2017_norm_bc-taxo.csv
- mOTU-v2: EC2M3-WGS-motus_rev-04-11-2019_bc-taxo.csv
- mOTU-v2: EC2M3-WGS-motus_rev-04-11-2019.csv

Data are normalized.
DIAMOND table is not used in the manuscript.

Taxonomy had sometime to be corrected because exports can generate, 
from some taxpath, empty fields which create a shift in the
taxonomical levels :

- script `detaxize2.py`
- script `taxo-bc-2.py`

Example :

```bash
$ python ./detaxize2.py > EC2M3-WGS-motus_rev-04-11-2019_bc-taxo.csv
$ python ./taxo-bc-2.py Comparison-EC2M3-WGS_Taxonomy-R1_rev-17-05-2017.csv > Comparison-EC2M3-WGS_Taxonomy-R1_rev-17-05-2017_bc-taxo.csv
```

Summarisation at genus level is done using the script `mOTU-summation-genus-keep-controls.R`. This
script keep only the 50 control individuals from the "French cohort" (mOTU-v2 was run on the
whole cohort).

