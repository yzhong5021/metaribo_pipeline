# MetaRibo-ITP processing and analysis pipeline

This repo contains the full preprocessing and analysis workflow for paired metagenomic and Ribo-ITP data.

## 1 - Environment setup

Build and activate the Conda environment:

```bash
conda env create -f envs/metaribo_pipeline_env.yaml
conda activate metaribo_pipeline
```

Apptainer is required by nf-core/mag, but is not in this Conda environment; confirm that Apptainer is on PATH before running metagenomic processing.

## 2 - Prepare config and sample files

Edit copies of the template configs stored in `config/` (`config/*_template.yaml):

```bash
cp config/master_template.yaml config/master.yaml
cp config/metagenome_preprocess_template.yaml config/metagenome_preprocess.yaml
cp config/contig_annotation_template.yaml config/contig_annotation.yaml
cp config/ribo_preprocess_template.yaml config/ribo_preprocess.yaml
cp config/samples_template.csv config/samples.csv
cp config/samples_ribo_template.csv config/samples_ribo.csv
```

- `config/master.yaml`: provide the repo root, scratch/result/log directories, file manifest, and the stage-specific config files.
- `file_manifest.json` is used to save and pass file paths between stages. Do not edit unless you are supplying precomputed files. Either provide an empty file or let the pipeline create the manifest itself.

### Specifying samples

- Samples should be in gzipped FASTQ format.
- `config/samples.csv`: metagenome sample sheet. For additional help, read the [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation. This path should be specified in `config/metagenome_preprocess.yaml`, under `samplesheet`.
    - Samples should contain the following entries: `sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform`
    - Sample ID's (`sample`) should be unique. Currently requires paired Illumina reads.
- `config/samples_ribo.csv`: Ribo-ITP sample sheet. This path should be specified in `config/ribo_preprocess.yaml`, under `samplesheet`.
    - Samples should contain `id,meta_id,path`.
    - Sample ID's (`id`) should be unique. `meta_id` should match the respective metagenome ID specified above in (`sample`).

## 3. Prepare reference databases

Download and unpack the CAT 2024-12-12 database:

```bash
mkdir -p /path/to/CAT
cd /path/to/CAT
wget -c https://tbb.bio.uu.nl/tina/CAT_pack_prepare/20241212_CAT_nr_website.tar.gz
tar -xzf 20241212_CAT_nr_website.tar.gz
```

Set `config/contig_annotation.yaml` paths:

```yaml
cat_db: /path/to/CAT/20241212_CAT_nr_website/db
cat_tax: /path/to/CAT/20241212_CAT_nr_website/tax
```

Download Bakta database:

```bash
mkdir -p /path/to/bakta/database
bakta_db download --output /path/to/bakta/database --type full
```

Set `config/contig_annotation.yaml` paths:

```yaml
bakta_db: /path/to/bakta/database/db
```

Download Kofam profiles:

```bash
mkdir -p /path/to/kofam
cd /path/to/kofam
wget -c https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget -c https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
gunzip ko_list.gz
tar -xzf profiles.tar.gz
```

Set `config/contig_annotation.yaml` paths:

```yaml
kofam_profiles: /path/to/kofam/profiles
kofam_list: /path/to/kofam/ko_list
```

Functional analyses also require a mapping of KO to KEGG terms (our study used BRITE-C); for example, `/path/to/ko_to_brite.tsv`. It should follow the below format:

```text
ko:K00001	Glycolysis / Gluconeogenesis
ko:K00001	Pyruvate metabolism
```

## 4. Run the pipeline

Run the three stages from the root. Logs for each stage are written to `{logs_dir}/metagenome_preprocess.log`, `{logs_dir}/contig_annotation.log`, and `{logs_dir}/ribo_preprocess.log`.

### Metagenome preprocessing and MAG generation:

- `config/params.mag.yaml`: nf-core/mag parameters. visit [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation for additional help.
- `config/nextflow.config`: nf-core/mag config. visit [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation for additional help.

```bash
./src/scripts/metagenome_preprocess.sh --master-config config/master.yaml
```

This runs nf-core/mag 5.4.2 with Apptainer and `-resume` (fastp, MEGAHIT, assembly backmapping via Bowtie2, six binners (MetaBAT2, MaxBin2, CONCOCT, MetaBinner, SemiBin2, and COMEBin), DAS Tool refinement, BUSCO and CheckM2, and GTDB-Tk classification)

The complete nf-core/mag results are output in `{scratch}/metagenome_preprocess`; the Nextflow work directory is stored in `{scratch}/work/metagenome_preprocess`. The relevant output and intermediate files are in `{results_dir}/metagenome`:

- Assembled contigs: `contigs/MEGAHIT-{META_ID}.contigs.fa`
- Bins: `bins/*.fa`
- Backmapping: `bowtie2/*.bam` and BAM indexes
- Contig taxonomy mapping: `contig_to_genus/contig_to_genus_{META_ID}.tsv`

### Contig annotation

```bash
./src/scripts/contig_annotation.sh --master-config config/master.yaml
```

This runs Bakta annotation, computes metagenome coverage for each CDS, runs KofamScan annotation, classifies remaining contigs via CAT, and initializes the master sheets for each metagenome sample.

Large intermediates are stored in `{scratch}/metagenome_annotation`; work files are in `{scratch}/work/metagenome_annotation`. The relevant output and intermediate files are in `{results_dir}/metagenome`:

- Bakta annotations: `bakta_cds_annotations_gff/{META_ID}_cds.gff3`
- ncRNA regions: `bakta_ncrna_bed/{META_ID}_ncrna.bed`
- CDS regions: `bakta_cds_bed/{META_ID}_cds.bed`
- CDS protein sequences: `bakta_cds_faa/{META_ID}_cds.faa`
- Metagenome coverage: `bakta_gene_mean_coverage/{META_ID}_gene_mean_coverage.tsv`
- KEGG annotations: `kofam_annotations/{META_ID}_kofam.tsv`
- Updated contig taxonomy mappings: `contig_to_genus/contig_to_genus_{META_ID}.tsv`
- Initialized master sheets (used downstream for analyses): `master_sheets/{META_ID}_master_sheet.tsv`

### Ribo-ITP preprocessing and alignment

```bash
./src/scripts/ribo_preprocess.sh --master-config config/master.yaml
```

This trims adapters and quality-filters reads, aligns to the metagenome contigs via Bowtie2, removes multimappers and ncRNA overlaps, and performs UMI deduplication.

Large intermediates are under `{scratch}/ribo_preprocess`. The relevant output and intermediate files are in `{results_dir}/ribo`:

- Final alignment file: `alignment/{RIBO_ID}.nomulti.filtered.sorted.deduped.bam` (and index)
- Alignment statistics: `stats/{RIBO_ID}.alignment_stats.txt`

## 5. Analyze outputs

The default notebooks assume the following path for the results directory: `results_dir: output`. If a different directory is intended, edit the paths in the I/O cell accordingly.

Set up and start JupyterLab:
```bash
python -m ipykernel install --user --name metaribo_pipeline --display-name metaribo_pipeline
jupyter lab
```

First, run `notebooks/coverage.ipynb`, which performs RPF binning and metagene calculation.

Then run `notebooks/functional_analysis.ipynb`, which performs read normalization, taxonomy/function profiling and visualization, and null permutation analysis.
