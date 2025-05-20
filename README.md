# RNA-seq pipeline

This repository contains an RNA-seq data processing pipeline implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). The pipeline automates the essential steps of an RNA-seq analysis, including Quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), Alignment ([STAR](https://github.com/alexdobin/STAR)), Quantification ([FeatureCount](https://subread.sourceforge.net/featureCounts.html)), and BedGraph generation.

## Prerequisites
Before running the pipeline, ensure that you have the following installed:
  - Conda or Miniconda
  - Snakemake v8.14.0

## Getting started
### 1. Clone the repository
  ```bash
  git clone https://github.com/krupas26/RNA-seq_pipeline.git
  cd RNA-seq_pipeline
  ```

### 2. Create a Conda Environment
   Use the provided env_rnaseq.yml to create the environment with all the required dependencies:
   ```
   conda env create --file env_rnaseq.yml
   ```
  Once the environment is created, activate it:
  ```
  conda activate rnaseq
  ```

### 3. Configure the Pipeline
   To keep sensitive data such as the file paths private, this project uses a configuration file that is excluded from version control.
   You can create a `config.yml` file in your favorite text editor and update the placeholder values with the paths and settings specific to your environment.

   #### Example snippet of `config.yml` file
   ```yaml
    workflow-profile: "/path/to/your/workflow/profile"
    metadata: "metadata.xlsx"
    chrom_sizes: "/path/to/chrom_sizes"
    adapter: "/path/to/adapter/file"
    fastq_dir: "/path/to/fastq/files"
    sample_column: "SampleID"
    samples_list: null
    run_id: "<your experiement identifier>"
  ```

### 4. Running the pipeline
  Run the pipeline with Snakemake:
  ```bash
  snakemake -s rnaseq_pipeline.smk --configfile <config.yml> --cores <number-of-cores>
  ```
##### To submit jobs on cluster as sbatch jobs
  Create a generic profile_config.yml
  ```yaml
    executor: "cluster-generic"
    cluster-generic-submit-cmd: "sbatch --output=logs/{rule}_{wildcards}_%j.out --error=logs/{rule}_{wildcards}_%j.err --mem={resources.mem_mb}"
    jobname: "{wildcards}_{rule}_{jobid}"
    default-resources: [mem_mb=4000, time='02:00:00', threads=1]
    use-conda: true
    jobs: 100
    rerun-incomplete: true
    latency-wait: 60
  ```
  Then run the pipeline with Snakemake as:
  ```bash
  snakemake -s rnaseq_pipeline.smk --profile <profile_config.yml> --configfile <config.yml> --cores <number-of-cores>
  ```




     
