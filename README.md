# MetaRibo-ITP Processing / Analysis Pipeline

This repo contains the full preprocessing and analysis workflow for paired metagenomic and Ribo-ITP data.

## 1 - Environment setup

Build and activate the single conda environment:

```bash
conda env create -f envs/metaribo_pipeline_env.yaml
conda activate metaribo_pipeline
```

## 2 - Prepare configs and samples

Use the templates provided in `config` as a guide.


`config/master.yaml`

- `root`: the path to this repo
- `scratch`: store large outputs here. Work dirs are always under `scratch/work`.
- `results_dir`: local results directory; can be relative to `root`, or the absolute path
- `logs_dir`: log directory; can be relative to `root`, or the absolute path
- `file_manifest`: the manifest the scripts update as each stage finishes
- `ribo_preprocess_params`, `metagenome_preprocess_params`, `contig_annotation_params`: stage-specific yaml paths; can be relative to `root`, or the absolute path

### Inputting Samples
- samples should be in gzipped FASTQ format.
- `config/samples.csv`: metagenome sample sheet. For help, visit the [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation. This path should be specified in `config/metagenome_preprocess.yaml`, under `samplesheet`.
    - Samples should contain `sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform,long_reads_platform`
- `config/samples_ribo.csv`: Ribo-ITP sample sheet. This path should be specified in `config/ribo_preprocess.yaml`, under `samplesheet`.
    - Samples should contain `id,meta_id,path`. `meta_id` should match the respective metagenome ID specified above in (`sample`).
- `config/metagenome_preprocess.yaml`: metagenome preprocessing-specific parameters.
- `config/params.mag.yaml`: nf-core/mag parameters. visit [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation for additional help.
- `config/nextflow.config`: nf-core/mag config. visit [nf-core/mag](https://nf-co.re/mag/5.4.2) documentation for additional help.
- `config/contig_annotation.yaml`: contig annotation specific parameters.
- `config/ribo_preprocess.yaml`: Ribo-ITP preprocessing/alignment parameters.

The intermediate files, outputs, and logs for specific stages are located in: 
- `{scratch}/{metagenome_preprocess,metagenome_annotation,ribo_preprocess}`; - `scratch/work/{metagenome_preprocess,metagenome_annotation,ribo_preprocess}`,  - `results_dir/{metagenome,ribo,coverage}`, and - `logs_dir/{metagenome_preprocess.log,contig_annotation.log,ribo_preprocess.log}`.

`file_manifest.json` is used to save and pass file paths between stages. Do not edit unless you are supplying precomputed files.

## 2 - Database generation

Download CAT database:

```bash
mkdir -p /path/to/cat/database
cd /path/to/cat/database
wget -c https://tbb.bio.uu.nl/tina/CAT_pack_prepare/20241212_CAT_nr_website.tar.gz
tar -xzf 20241212_CAT_nr_website.tar.gz
```

Set `config/contig_annotation.yaml` paths:

```yaml
cat_db: /path/to/cat/database/20241212_CAT_nr_website/db
cat_tax: /path/to/cat/database/20241212_CAT_nr_website/tax
```

Download Bakta database:

```bash
mkdir -p /path/to/bakta/database
bakta_db download --output /path/to/bakta/database --type full
```

Set `config/contig_annotation.yaml` paths:

```yaml
bakta_db: /path/to/bakta/database/
```

Download Kofam profiles:

```bash
mkdir -p /path/to/kofam/database/
cd /path/to/kofam/database/
wget -c https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget -c https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
gunzip ko_list.gz
tar -xzf profiles.tar.gz
```

Set `config/contig_annotation.yaml` paths:

```yaml
kofam_profiles: /path/to/kofam/database/profiles
kofam_list: /path/to/kofam/database/ko_list
```

## 3 - Run the pipeline

Run each stage from the root.

Run metagenome preprocessing / MAG generation:

```bash
./src/scripts/metagenome_preprocess.sh --master-config config/master.yaml
```

This runs nf-core/mag, performs backmapping and contig genus assignment, then collects relevant outputs (contigs, BAMs, bins, and contig genus assignments) into the specified local output.

Run contig annotation:

```bash
./src/scripts/contig_annotation.sh --master-config config/master.yaml
```

This runs Bakta, CAT, metagenome coverage, KofamScan, and initializes sample-specific CDS statistics sheets.

Run Ribo-ITP preprocessing and alignment:

```bash
./src/scripts/ribo_preprocess.sh --master-config config/master.yaml
```

This trims reads, extracts UMIs when enabled, builds the Bowtie2 index from the matched metagenome contigs, aligns ribo reads, removes multimappers and ncRNA overlaps when enabled, deduplicates UMIs when enabled, writes alignment stats, and updates `file_manifest.json`.

## 4 - Analyze outputs

Set up JupyterLab:

```bash
python -m ipykernel install --user --name=metaribo_pipeline --display-name "metaribo_pipeline"
jupyter lab
```

Run `notebooks/coverage.ipynb` first: ribo read binning and metagene calculation;

Then run `notebooks/functional_analysis.ipynb`: read normalization, taxonomy/functional visualization, null permutation analysis.

## Outputs

Logs are written under `logs_dir`.

Local outputs are written under `results_dir`. Large outputs are written under `scratch`; work directories use `scratch/work`.

Essential intermediate files are stored in `file_manifest.json`.
