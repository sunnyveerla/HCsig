__authors__ = "Srinivas Veerla, Guyuan TANG"
__copyright__ = "Copyright 2024, Srinivas Veerla"
__email__ = "srinivas.veerla@med.lu.se"
__license__ = "GPL-3"

import pandas as pd

#Configuration file
configfile: "config/config.yaml"

#Read sample information file
sample_info = (pd.read_csv(config['samples'],     
    dtype={'sample_name':str, 'patient':str, 'type':str,'fastq_1':str, 'fastq_2':str})
    .set_index('sample_name', drop=False))

#Working directory location
work_dir = config['workdir']

#Result data files storage location
results = config['outputdir']

rule all:	
    input:  expand(results + '02_alignment/{sample}/{sample}.sorted.bam.bai', sample=sample_info.sample_name)	

########### 1 sWGS raw data preprocessing ##########
rule fastp:
    input:
        R1 = lambda wildcards: sample_info.loc[wildcards.sample, 'fastq_1'],
        R2 = lambda wildcards: sample_info.loc[wildcards.sample, 'fastq_2']
    output:
        R1 = results + '01_preprocess/{sample}/reads/{sample}_R1_preprocess.fastq.gz',
        html = results + '01_preprocess/{sample}/html/{sample}_fastp.html',
        R2 = results + '01_preprocess/{sample}/reads/{sample}_R2_preprocess.fastq.gz' 
    log: 'log/fastp/{sample}_fastp.log'
    threads: 20
    params: json = results + '01_preprocess/{sample}/html/{sample}_fastp.json'
    #conda: "envs/preprocess_env.yaml"
    shell: """
        fastp --detect_adapter_for_pe \
        --correction --cut_right --thread {threads} \
        --html {output.html} --json {params.json} \
        --in1 {input.R1} --in2 {input.R2} \
        --out1 {output.R1} --out2 {output.R2} \
        2>{log}
    
    rm {params.json}
    """
 # 1.1 quality assessment of preprocessed reads with fastqc
rule fastqc:
    ### we will use fastqc to generate the quality control stats from the outputs of fastp
    input:
        R1_seq = results + '01_preprocess/{sample}/reads/{sample}_R1_preprocess.fastq.gz',
        R2_seq = results + '01_preprocess/{sample}/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        R1_html = results + '01_preprocess/{sample}/html/{sample}_R1_preprocess_fastqc.html',
        R1_qc = results + '01_preprocess/{sample}/reports/{sample}_R1_preprocess_fastqc.zip',
        R2_html = results + '01_preprocess/{sample}/html/{sample}_R2_preprocess_fastqc.html',
        R2_qc = results + '01_preprocess/{sample}/reports/{sample}_R2_preprocess_fastqc.zip'
    log: 'log/fastqc/{sample}.fastqc.log'
    params: 
        outdir = results + '01_preprocess/{sample}/reports/',
        out_html = results + '01_preprocess/{sample}/html/'
    threads: 20
    #conda: 'envs/preprocess_env.yaml'
    shell: """
    fastqc -o {params.outdir} {input.R1_seq} {input.R2_seq} 2>{log}
    mv {params.outdir}*_fastqc.html {params.out_html}
    """

########## 2 Alignment ####################

# 2 mapping the reads to the indexed reference genome
rule map_reads:
    ### use bwa again for alignment
    input: 
        #idx = rules.bwa_index.output,
        link_up = rules.fastqc.input,
        R1 = results + '01_preprocess/{sample}/reads/{sample}_R1_preprocess.fastq.gz',
        R2 = results + '01_preprocess/{sample}/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        results + '02_alignment/{sample}/{sample}.unsorted.sam'
    log: 'log/bwa_mapping/{sample}/{sample}.log'
    params:
        index_ref = 'resources/genome/hg38.fa'
    #conda: 'envs/alignment.yaml'
    threads: config['bwa_mapping']['threads']
    shell: """
    bwa mem -M -t {threads} \
        {params.index_ref} {input.R1} {input.R2} > {output} \
        2>{log}
    """

# 2.1 sorting the SAM files
rule sort_sam: 
    ### using Picard to sort the sam files 
    input:
        results + '02_alignment/{sample}/{sample}.unsorted.sam'        
    output:
        results + '02_alignment/{sample}/{sample}.sorted.sam'
    log: 'log/sort_sam/{sample}.log'
    #conda: 'envs/clean_up.yaml'
    shell: """
    picard SortSam \
        INPUT={input} \
        OUTPUT={output} \
        SORT_ORDER=coordinate \
        2>{log}
    """

#2.2 Convert SAM  file to BAM file 
rule convert_sam_2_bam:
    input:
        results + '02_alignment/{sample}/{sample}.sorted.sam'
    output:
        results + '02_alignment/{sample}/{sample}.sorted.bam'
    log: 'log/sam_2_bam/{sample}.log'
    threads:20
    shell: """
    samtools view -@ {threads} -bo {output} {input}
    """ 

#2.3  indexing the BAM files and removing unsorted and sorted sam files
rule index_bam:
    ### using samtools to show the stats of the sorted and deduplicates outputs and to index the bam files
    input:
       bam = results + '02_alignment/{sample}/{sample}.sorted.bam',
       del_unsort_sam = results + '02_alignment/{sample}/{sample}.unsorted.sam',
       del_sort_sam = results + '02_alignment/{sample}/{sample}.sorted.sam'
    output:
       results + '02_alignment/{sample}/{sample}.sorted.bam.bai'    
    log: 'log/bam_stat/{sample}.log'
    threads: 20
    #conda: 'envs/clean_up.yaml'
    shell: """
     samtools flagstat {input.bam} | tee {log}    
     samtools index -@ {threads} {input.bam} 
     rm -f {input.del_unsort_sam}
     rm -f {input.del_sort_sam}    
    #qualimap bamqc -bam {input} --java-mem-size=4G
    """


