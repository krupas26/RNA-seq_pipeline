#import necessary libraries
import pandas as pd
import os
import subprocess
import yaml
import re

#load variables from config file
filepath = config['workflow-profile']
print("Filepath: ", filepath)
metadata = config['metadata_path']
chrom_sizes = config['chrom_sizes_path']
adapter = config['adapter_file_path']
fastq_dir = config.get('fastq_dir', os.path.join(filepath, 'data/fastq'))
refGenome_dir = config['refGenome_dir_path']
gtf_dir = config['gtf_dir_path']

#either parse sample column from metadata or use the provided sample list directly
sample_column = config.get('sample_column', None)
samples_list = config.get('samples_list', None)

#extract samples
if samples_list:
	print('Using provided sample list for the analysis..')
	samples = samples_list
else:
	if not metadata or not sample_column:
		raise ValueError("Metadata file or sample column not provided. Please add `metadata` to config file or add `sample_column` to the config file.")
	
	#define path for metadata
	file = os.path.join(filepath, metadata)

	#extract the sample names from file and make sample list
	df = pd.read_excel(file)
	if sample_column not in df.columns:
		raise ValueError(f"Column '{sample_column}' not found in the metadata file.")

	samples = df[sample_column].tolist()

print(f"Samples: {samples}")

#experiment ID - used for naming the files
run_id = config.get('run_id', None)
print("Run id : ", run_id)

rule all:
	input:
		expand(f"{filepath}/qc/{{sample}}/fastqc/{{sample}}_{{read}}_fastqc{{ext}}", sample=samples, read=[1,2], ext=[".html", ".zip"]),
		expand(f"{filepath}/qc/{{sample}}/fastq_screen/{{sample}}_{{read}}_screen.html", sample=samples, read=[1,2]),
		expand(f"{filepath}/data/trimmed_fastq/{{sample}}_{{read}}.fastq.gz", sample=samples, read=['1P', '2P', '1U', '2U']),
		expand(f"{filepath}/qc/{{sample}}/trimmed_fastqc/{{sample}}_{{read}}_fastqc.html", sample=samples, read=['1P', '2P']),
		expand(f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=samples),
		expand(f"{filepath}/qc/picard_metrics/{{sample}}_marked_dups.txt", sample=samples),
		expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups.bam", sample=samples),
		expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups_flagstats.txt", sample=samples),
		expand(f"{filepath}/data/aligned/{{sample}}/{{sample}}_sorted.bam.bai", sample=samples),
		f"{filepath}/{run_id}_bam_files.txt",
		f"{filepath}/{run_id}_counts.txt",
		directory(f"{filepath}/{run_id}_multiqc_out"),
		expand(f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.out.bg", sample=samples),
		expand(f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.out.bg", sample=samples),
		expand(f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.signFlipped.sorted.out.bg", sample=samples),
		expand(f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.sorted.out.bg", sample=samples),
		expand(f"{filepath}/data/tracks/{{sample}}_str1_minus.bw", sample=samples),
		expand(f"{filepath}/data/tracks/{{sample}}_str2_plus.bw", sample=samples)
	localrule: True

#Data Preprocessing
#Quality control using FastQC and Fastq Screen
rule fastqc:
	input:
		fq1 = f"{fastq_dir}/{{sample}}_1.fastq.gz",
		fq2 = f"{fastq_dir}/{{sample}}_2.fastq.gz"
	output:
		multiext(f"{filepath}/qc/{{sample}}/fastqc/{{sample}}_1_fastqc", ".html", ".zip"),
		multiext(f"{filepath}/qc/{{sample}}/fastqc/{{sample}}_2_fastqc", ".html", ".zip")
	params:
		dir= f"{filepath}/qc/{{sample}}/fastqc/"
	conda:
		"rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		fastqc {input.fq1} {input.fq2} -o {params.dir}
		"""

rule fastq_screen:
	input:
		fq1 = f"{filepath}/data/fastq/{{sample}}_1.fastq.gz",
		fq2 = f"{filepath}/data/fastq/{{sample}}_2.fastq.gz"
	output:
		multiext(f"{filepath}/qc/{{sample}}/fastq_screen/{{sample}}_1_screen", ".html", ".png"),
		multiext(f"{filepath}/qc/{{sample}}/fastq_screen/{{sample}}_2_screen", ".html", ".png")
	params:
		outdir = f"{filepath}/qc/{{sample}}/fastq_screen/"
	conda:
		"rnaseq"
	resources:
		threads=32,
		time="12:00:00",
		mem_mb=100000
	shell:
		"""
		fastq_screen --aligner bwa --outdir {params.outdir} {input.fq1} {input.fq2} --force --threads {resources.threads}
		"""
		
#Trim the fastq sequences and removing adapter contamination using Trimmomatic ---- used Trim_galore (if trimmomatic fails to work). Follow the script run_trim_galore.sh
rule trimmomatic:
	input:
		fq1 = f"{filepath}/data/fastq/{{sample}}_1.fastq.gz",
		fq2 = f"{filepath}/data/fastq/{{sample}}_2.fastq.gz"
	output:
		trimmed_fq1 = f"{filepath}/data/trimmed_fastq/{{sample}}_1P.fastq.gz",
		trimmed_fq2 = f"{filepath}/data/trimmed_fastq/{{sample}}_2P.fastq.gz",
		unpaired_fq1 = f"{filepath}/data/trimmed_fastq/{{sample}}_1U.fastq.gz",
		unpaired_fq2 = f"{filepath}/data/trimmed_fastq/{{sample}}_2U.fastq.gz"
	params:
		outdir = f"{filepath}/data/trimmed_fastq/",
		adapter = f"{adapter}"
	conda:
		"rnaseq"
	resources:
		threads=32,
		mem_mb=50000,
		time="24:00:00"
	shell:
		"""
		trimmomatic PE -threads {resources.threads} -phred33 \
		{input.fq1} {input.fq2} {output.trimmed_fq1} {output.unpaired_fq1} {output.trimmed_fq2} {output.unpaired_fq2} \
		ILLUMINACLIP:{params.adapter}:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:10 MINLEN:36
		"""
		
#Repeat FastQC to confirm the quality ----step performed using Trim_galore (if Trimmomatic fails 
rule re_fastqc:
	input:
		trimmed_fq1 = f"{filepath}/data/trimmed_fastq/{{sample}}_1P.fastq.gz",
		trimmed_fq2 = f"{filepath}/data/trimmed_fastq/{{sample}}_2P.fastq.gz"
	output:
		multiext(f"{filepath}/qc/{{sample}}/trimmed_fastqc/{{sample}}_1P_fastqc", ".html", ".zip"),
		multiext(f"{filepath}/qc/{{sample}}/trimmed_fastqc/{{sample}}_2P_fastqc", ".html", ".zip")
	params:
		dir= f"{filepath}/qc/{{sample}}/trimmed_fastqc/"
	conda:
		"rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		fastqc {input.trimmed_fq1} {input.trimmed_fq2} -o {params.dir}
		"""
		
#Alignment using STAR
rule star:
	input:
		trimmed_fq1 = f"{filepath}/data/trimmed_fastq/{{sample}}_1P.fastq.gz",
		trimmed_fq2 = f"{filepath}/data/trimmed_fastq/{{sample}}_2P.fastq.gz"
	output:
		bam_files = f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
	params:
		ref = f"{refGenome_dir}",
		outdir = f"{filepath}/data/aligned/{{sample}}/"
	conda:
		"rnaseq"
	resources:
		threads=32,
		mem_mb=70000,
		time="24:00:00"
	shell:
		"""
		#increase the number of open files at a time
		ulimit -n 4096
		STAR --runThreadN {resources.threads} --genomeDir {params.ref} --readFilesCommand zcat \
		--readFilesIn {input.trimmed_fq1} {input.trimmed_fq2} \
		--outSAMtype BAM SortedByCoordinate --outFileNamePrefix "{params.outdir}{wildcards.sample}_"
		"""
		
#Mark duplicates using Picard MarkDuplicates
rule mark_duplicates:
	input:
		bam_files = f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		metrics = f"{filepath}/qc/picard_metrics/{{sample}}_marked_dups.txt",
		marked_dup_bam = f"{filepath}/data/picard_bam/{{sample}}_marked_dups.bam"
	conda: "rnaseq"
	resources:
		threads=16,
		time="12:00:00",
		mem_mb=160000
	shell:
		"""
		picard MarkDuplicates I={input.bam_files} O={output.marked_dup_bam} M={output.metrics} \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
		"""
		
#Get the statistics of alignment using samtools flagstats
rule flagstats:
	input:
		bam_files = f"{filepath}/data/picard_bam/{{sample}}_marked_dups.bam"
		#bam_files = f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		bam_stats = f"{filepath}/data/picard_bam/{{sample}}_marked_dups_flagstats.txt"
	conda: 
		"rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		samtools flagstat {input.bam_files} > {output.bam_stats}
		"""

#Index the bam files
rule index_bam:
	input:
		bam_files = f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		sorted_indexed_bam = f"{filepath}/data/aligned/{{sample}}/{{sample}}_sorted.bam.bai"
	conda: 
		"rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		samtools index {input.bam_files} {output.sorted_indexed_bam}
		"""

#Combine bam files to a text file
rule combine_bam:
	input:
		bam_files = expand(f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=samples)
	output:
		bam_txt = f"{filepath}/{run_id}_bam_files.txt"
	localrule: True
	shell:
		"""
		#find {filepath}/data/aligned/ -type f -name "*.bam" | sort > {output.bam_txt}
		for sample in {input.bam_files}; do
					echo $sample >> {output.bam_txt}
		done
		"""

#Extract counts using featureCounts
rule counts:
	input:
		bam_files = f"{filepath}/{run_id}_bam_files.txt",
		gtf_file = f"{gtf_dir}/GRCh38v29+ERCC.gtf"
	output:
		counts = f"{filepath}/{run_id}_counts.txt",
		counts_summary = f"{filepath}/{run_id}_counts.txt.summary"
	conda: 
		"rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		#xargs -a {input.bam_files} -I {{}} featureCounts -p --countReadPairs -T 32 -a {input.gtf_file} -o {output.counts} {{}}
		featureCounts -p --countReadPairs -s 2 -T 32 -a {input.gtf_file} -o {output.counts} $(cat {input.bam_files})
		"""

#Combine all reports using MultiQC
rule multiqc:
	input: expand(f"{filepath}/qc/{{sample}}/fastqc/{{sample}}_1_fastqc.html", sample=samples),
		expand(f"{filepath}/qc/{{sample}}/fastqc/{{sample}}_2_fastqc.html", sample=samples),
		expand(f"{filepath}/qc/{{sample}}/fastq_screen/{{sample}}_1_screen.html", sample=samples),
		expand(f"{filepath}/qc/{{sample}}/fastq_screen/{{sample}}_2_screen.html", sample=samples),
		expand(f"{filepath}/data/trimmed_fastq/{{sample}}_1P.fastq.gz", sample=samples),
		expand(f"{filepath}/data/trimmed_fastq/{{sample}}_2P.fastq.gz", sample=samples),
		expand(f"{filepath}/data/trimmed_fastq/{{sample}}_1U.fastq.gz", sample=samples),
		expand(f"{filepath}/data/trimmed_fastq/{{sample}}_2U.fastq.gz", sample=samples),
		expand(f"{filepath}/qc/{{sample}}/trimmed_fastqc/{{sample}}_1P_fastqc.html", sample=samples),
		expand(f"{filepath}/qc/{{sample}}/trimmed_fastqc/{{sample}}_2P_fastqc.html", sample=samples),
		expand(f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=samples),
		expand(f"{filepath}/qc/picard_metrics/{{sample}}_marked_dups.txt", sample=samples),
		expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups.bam", sample=samples),
		expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups_flagstats.txt", sample=samples),
		expand(f"{filepath}/data/aligned/{{sample}}/{{sample}}_sorted.bam.bai", sample=samples),
		f"{filepath}/{run_id}_counts.txt.summary"
	output: directory(f"{filepath}/{run_id}_multiqc_out")
	params:
		qc_dir = expand(f"{filepath}/qc/{{sample}}", sample=samples),
		aligned_dir = expand(f"{filepath}/data/aligned/{{sample}}", sample=samples),
		picard_metrics = expand(f"{filepath}/qc/picard_metrics/{{sample}}_marked_dups.txt", sample=samples),
		picard_bam = expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups.bam", sample=samples),
		picard_bam_stats = expand(f"{filepath}/data/picard_bam/{{sample}}_marked_dups_flagstats.txt", sample=samples),
		counts = f"{filepath}/{run_id}_counts.txt.summary"
	conda: 'rnaseq'
	resources: 
		threads=8,
		time="1:00:00"
	localrule: True
	shell:
		"""
		multiqc {params.qc_dir} {params.aligned_dir} {params.picard_metrics} {params.picard_bam} {params.picard_bam_stats} {params.counts} -o {output}
		"""

#Make Tracks -- uncomment when required (Don't forget to uncomment output files from rule all as well!!)
#Get bedGraph files using STAR
rule bedgraph:
	input: 
		bam_file = f"{filepath}/data/aligned/{{sample}}/{{sample}}_Aligned.sortedByCoord.out.bam"
	output:
		bg_str1 = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.out.bg",
		bg_str2 = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.out.bg"
	params:
		outdir = f"{filepath}/data/bedgraph"
	conda: "rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		STAR --runMode inputAlignmentsFromBAM --inputBAMfile {input.bam_file} \
		--outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM --outFileNamePrefix {params.outdir}/{wildcards.sample}_
		"""

#Sort the bedgraph
rule sort_bg:
	input:
		bg_str1 = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.out.bg",
		bg_str2 = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.out.bg"
	output:
		bg1_sorted = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.signFlipped.sorted.out.bg",
		bg2_sorted = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.sorted.out.bg"
	conda: "rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		awk -F '\\t' 'BEGIN {{OFS = FS}} {{$4 = $4 * -1; print $0}}' {input.bg_str1} | LC_COLLATE=C sort -S 7G -k 1,1 -k 2,2n > {output.bg1_sorted}

		LC_COLLATE=C sort -S 7G -k 1,1 -k 2,2n {input.bg_str2} >  {output.bg2_sorted}
		"""

#Generate the bigwig tracks
rule tracks:
	input:
		bg1_sorted = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str1.signFlipped.sorted.out.bg",
		bg2_sorted = f"{filepath}/data/bedgraph/{{sample}}_Signal.Unique.str2.sorted.out.bg"
	output:
		bw_str1 = f"{filepath}/data/tracks/{{sample}}_str1_minus.bw",
		bw_str2 = f"{filepath}/data/tracks/{{sample}}_str2_plus.bw"
	params:
		chrom_sizes = f"{chrom_sizes}"
	conda: "rnaseq"
	resources:
		threads=8,
		time="12:00:00",
		mem_mb=10000
	shell:
		"""
		bedGraphToBigWig {input.bg1_sorted} {params.chrom_sizes} {output.bw_str1}
		bedGraphToBigWig {input.bg2_sorted} {params.chrom_sizes} {output.bw_str2}
		"""
