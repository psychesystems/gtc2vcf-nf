# gtc2vcf-nf

Nextflow workflow for running [gtc2vcf](https://github.com/freeseek/gtc2vcf).

## Example run

### Identify chip

```sh
nextflow run workflows/gtc2vcf.nf -resume -qs 8 --command identify --idats "data/idats/*.idat"
```
```sh
Chip type guess (G/R): [GSA-MD-48v4-0_20098041, GSA-MD-48v4-0_20098041] Count: 576
```

Download manifest for chip:

```sh
mkdir data/manifest
wget -P data/manifest "https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/InfiniumGlobalScreeningArrayv4.0ClusterFile-egt.zip" 
unzip data/manifest/InfiniumGlobalScreeningArrayv4.0ClusterFile-egt.zip -d data/manifest
wget -P data/manifest https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/GSA-48v4-0_20085471_D2.bpm
wget -P data/manifest https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/GSA-48v4-0_20085471_D2.csv
```

### Convert to GTC

```sh
nextflow run workflows/gtc2vcf.nf -resume -qs 8 --command gtc \
--idats "data/idats/*.idat" \
--egt data/manifest/GSA-48v4-0_20085471_D2_ClusterFile.egt \
--bpm data/manifest/GSA-48v4-0_20085471_D2.bpm
```

GTC files are written to `data/gtc/GSA-48v4-0_20085471_D2`.

### Convert to VCF

Download and uncompress FASTA reference files.

```sh
mkdir data/fasta
wget -P data/fasta https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip data/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Index the FASTA and convert GTCs to VCF
```sh
nextflow run workflows/gtc2vcf.nf -resume -qs 8 --command vcf \
--gtc "data/gtc/GSA-48v4-0_20085471_D2/*.gtc" \
--egt data/manifest/GSA-48v4-0_20085471_D2_ClusterFile.egt \
--bpm data/manifest/GSA-48v4-0_20085471_D2.bpm \
--csv data/manifest/GSA-48v4-0_20085471_D2.csv \
--fasta data/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
--