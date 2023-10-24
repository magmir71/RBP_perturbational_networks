import os
import yaml
from pathlib import Path
import pandas as pd
from itertools import product
from functools import reduce 
import operator

localrules: all, prefetch_se, prefetch_pe

### SAMPLES TABLE ###
samples_table = pd.read_csv(config['samples'], header=0, index_col=False, sep=",")

samples_table.loc[samples_table['LibraryLayout'] == 'PAIRED', 'LibraryLayout'] = 'Paired'
samples_table.loc[samples_table['LibraryLayout'] == 'SINGLE', 'LibraryLayout'] = 'Single'

samples_table["absPath"] = config['outdir'] + os.path.sep + samples_table['LibraryLayout'] + os.path.sep + samples_table["sample"]
samples_table["absPathPrefix"] = samples_table["absPath"] + os.path.sep + samples_table["sample"]

END = ["_1","_2"]
ext = [".fastq.gz"]
short_ext = [".fastq"]

### PAIRED-END DICT PREP ###
PE_table = samples_table[samples_table["LibraryLayout"] == "Paired"]

Run_End = [run+end for run, end in product(PE_table["sample"], END)]

PE_path_fastq = [prefix+end+fastq for prefix, end, fastq in product(PE_table["absPathPrefix"], END, short_ext)]
PE_paths_gz = [prefix+end+gz for prefix, end, gz in product(PE_table["absPathPrefix"], END, ext)]

PE_fastq = pd.DataFrame({"Run_end":Run_End, "fastq":PE_path_fastq, "Path":PE_paths_gz, "Run":[s.split("_")[0] for s in Run_End]})

PE_nested_gz = {
    run: g_df.groupby('Run_end')["Path"].agg(list).to_dict()
    for run, g_df in PE_fastq.groupby("Run")
}

PE_nested_fastq = {
    run: g_df.groupby('Run_end')["fastq"].agg(list).to_dict()
    for run, g_df in PE_fastq.groupby("Run")
}


### SINGLE-END TABLE PREP ###
SE_table = samples_table[samples_table["LibraryLayout"] == "Single"]
SE_paths = [prefix+ext for prefix, ext in product(SE_table["absPathPrefix"], ext)]

SE_fastq = pd.DataFrame({"Run":SE_table["sample"], "Path":SE_paths})


####################
##### RULE ALL #####
####################

rule all:
    input:
        prefetch_se = expand(os.path.join(
            config["outdir"],
            "Single",
             "{run}", 
             "{run}.sra"),
             run = list(SE_fastq["Run"])),

        prefetch_pe = expand(os.path.join(
            config["outdir"],
            "Paired",
             "{run}", 
             "{run}.sra"),
             run = list(PE_fastq["Run"])),

        fasterq_dump_se = expand(os.path.join(
            config["outdir"], 
            "Single", 
            "{run}", 
            "{run}.dumped"),
             run = list(SE_fastq["Run"])),

        fasterq_dump_pe = expand(os.path.join(
            config["outdir"], 
            "Paired", 
            "{run}", 
            "{run}.dumped"),
            run = list(PE_fastq["Run"])),

        compress_fastq_se = expand(os.path.join(
            config['outdir'],
            "Paired",
             "{run}", 
             "{run}.fastq.gz"),
             run = list(SE_fastq['Run'])),

        compress_fastq_pe = expand(os.path.join(
            config['outdir'],
            "Paired",
            "{run}", 
            "{run}{end}.fastq.gz"),
            run = list(PE_fastq["Run"]),
            end = END),

        seqtk_se = expand(os.path.join(
            config['outdir'], 
            "Single",
            "{run}",
            "seqtk",
            "{run}.qc"),
            run = list(SE_fastq["Run"])),
        
        seqtk_pe = expand(os.path.join(
            config['outdir'], 
            "Paired", 
            "{run}", 
            "seqtk", 
            "{run}{end}.qc"),
            run = list(PE_fastq["Run"]),
            end = END),
               
        fastqc_se = expand(os.path.join(
            config['outdir'], 
            "Single",
            "{run}",
            "fastqc"),
            run = list(SE_fastq["Run"])),
        
        fastqc_pe = expand(os.path.join(
            config['outdir'],
            "Paired",
            "{run}",
            "fastqc",
            "{run}{end}"),
            run = list(PE_fastq["Run"]),
            end = END),
        
        star_nogtf_se = expand(os.path.join(
            config['outdir'],
            "Single",
            "{run}",
            "BAMS",
            "noGTF",  
            "{run}.Aligned.sortedByCoord.out.bam"),
            run = list(SE_fastq["Run"])),
        
        star_gtf_se = expand(os.path.join(
            config['outdir'],
            "Single",
            "{run}",
            "BAMS",
            "withGTF", 
            "{run}.Aligned.sortedByCoord.out.bam"),
            run = list(SE_fastq["Run"])),
                  
        star_nogtf_pe = expand(os.path.join(
            config['outdir'],
            "Paired",
            "{run}",
            "BAMS",
            "noGTF", 
            "{run}.Aligned.sortedByCoord.out.bam"),
            run = list(PE_fastq["Run"])),
        
        star_gtf_pe = expand(os.path.join(
            config['outdir'],
            "Paired",
            "{run}",
            "BAMS",
            "withGTF",  
            "{run}.Aligned.sortedByCoord.out.bam"),
            run = list(PE_fastq["Run"])),
        
        samtools_index_nogtf_se = expand(os.path.join(
            config['outdir'],
            "Single", 
            "{run}",
            "BAMS",
            "noGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.bai"),
            run = list(SE_fastq["Run"])),

        samtools_idxstats_nogtf_se = expand(os.path.join(
            config['outdir'],
            "Single",  
            "{run}",
            "BAMS",
            "noGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.idxstats"),
            run = list(SE_fastq["Run"])),
        
        samtools_index_gtf_se = expand(os.path.join(
            config['outdir'],
            "Single", 
            "{run}",
            "BAMS",
            "withGTF",
            "samtools",  
            "{run}.Aligned.sortedByCoord.out.bai"),
            run = list(SE_fastq["Run"])),

        samtools_idxstats_gtf_se = expand(os.path.join(
            config['outdir'], 
            "Single", 
            "{run}",
            "BAMS",
            "withGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.idxstats"),
            run = list(SE_fastq["Run"])),
        
        samtools_index_nogtf_pe = expand(os.path.join(
            config['outdir'], 
            "Paired",
            "{run}",
            "BAMS",
            "noGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.bai"),
            run = list(PE_fastq["Run"])),

        samtools_idxstats_nogtf_pe = expand(os.path.join(
            config['outdir'], 
            "Paired",
            "{run}",
            "BAMS",
            "noGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.idxstats"),
            run = list(PE_fastq["Run"])),
        
        samtools_index_gtf_pe = expand(os.path.join(
            config['outdir'],
            "Paired",
            "{run}",
            "BAMS",
            "withGTF",
            "samtools",  
            "{run}.Aligned.sortedByCoord.out.bai"),
            run = list(PE_fastq["Run"])),

        samtools_idxstats_gtf_pe = expand(os.path.join(
            config['outdir'],
            "Paired", 
            "{run}",
            "BAMS",
            "withGTF",
            "samtools", 
            "{run}.Aligned.sortedByCoord.out.idxstats"),
            run = list(PE_fastq["Run"])),

        megadepth_gtf_se_genes = expand(os.path.join(
            config['outdir'], 
            "Single", 
            "{run}", 
            "BAMS", 
            "withGTF", 
            "megadepth",
            "genes", 
            "{run}_auc.out"),
            run = list(SE_fastq["Run"])),
        
        megadepth_nogtf_se_genes = expand(os.path.join(
            config['outdir'], 
            "Single", 
            "{run}", 
            "BAMS", 
            "noGTF", 
            "megadepth",
            "genes", 
            "{run}_auc.out"),
            run = list(SE_fastq["Run"])),
        
        megadepth_gtf_se_exons = expand(os.path.join(
            config['outdir'], 
            "Single", 
            "{run}", 
            "BAMS", 
            "withGTF", 
            "megadepth",
            "exons", 
            "{run}_auc.out"),
            run = list(SE_fastq["Run"])),
        
        megadepth_nogtf_se_exons = expand(os.path.join(
            config['outdir'], 
            "Single", 
            "{run}", 
            "BAMS", 
            "noGTF", 
            "megadepth",
            "exons", 
            "{run}_auc.out"),
            run = list(SE_fastq["Run"])),

        megadepth_gtf_pe_genes = expand(os.path.join(
            config['outdir'], 
            "Paired", 
            "{run}", 
            "BAMS", 
            "withGTF", 
            "megadepth",
            "genes", 
            "{run}_auc.out"),
            run = list(PE_fastq["Run"])),
        
        megadepth_nogtf_pe_genes = expand(os.path.join(
            config['outdir'], 
            "Paired", 
            "{run}", 
            "BAMS", 
            "noGTF", 
            "megadepth",
            "genes", 
            "{run}_auc.out"),
            run = list(PE_fastq["Run"])),
        
        megadepth_gtf_pe_exons = expand(os.path.join(
            config['outdir'], 
            "Paired", 
            "{run}", 
            "BAMS", 
            "withGTF", 
            "megadepth",
            "exons", 
            "{run}_auc.out"),
            run = list(PE_fastq["Run"])),
        
        megadepth_nogtf_pe_exons = expand(os.path.join(
            config['outdir'], 
            "Paired", 
            "{run}", 
            "BAMS", 
            "noGTF", 
            "megadepth",
            "exons", 
            "{run}_auc.out"),
            run = list(PE_fastq["Run"])),

        DGEA = os.path.join(config['outdir'], "DEA", "prova.flag")

## custom functions ##
def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)

def getFastq1FromRunList(dict, run):
    return dict.get(run, {}).get(run+'_1')

def getFastq2FromRunList(dict, run):
    return dict.get(run, {}).get(run+'_2')



#### SRA PREFETCH - FASTQ DUMP - ZIP ####
rule prefetch_se:
    output:
        os.path.join(config['outdir'], "Single", "{run}", "{run}.sra"),
    params:
        prefetch = config['prefetch'],
        outdir = os.path.join(config['outdir'], "Single"),
    conda:
        "/scicore/home/zavolan/franch0002/myscripts/sra_tools.yaml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "prefetch.log")
    shell:
        """
        mkdir -p {params.outdir}; {params.prefetch} {wildcards.run} --output-directory {params.outdir} --force ALL &> {log}
        """

rule prefetch_pe:
    output:
        os.path.join(config['outdir'], "Paired", "{run}", "{run}.sra"),
    params:
        prefetch = config['prefetch'],
        outdir = os.path.join(config['outdir'], "Paired"),
    conda:
        "/scicore/home/zavolan/franch0002/myscripts/sra_tools.yaml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "prefetch.log")
    shell:
        """
        mkdir -p {params.outdir}; {params.prefetch} {wildcards.run} --output-directory {params.outdir} --force ALL &> {log}
        """

rule fasterq_dump_se:
    input:
        os.path.join(config['outdir'], "Single", "{run}", "{run}.sra"),
    output:
        flag = os.path.join(config['outdir'], "Single", "{run}", "{run}.dumped"),
        fq = os.path.join(config['outdir'], "Single", "{run}", "{run}.fastq"),
    params:
        outdir=os.path.join(config['outdir'], "Single", "{run}"),
        fastqdump = config['fasterq-dump'],
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    threads: 4
    conda:
        "/scicore/home/zavolan/franch0002/myscripts/sra_tools.yaml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "fasterq_dump.log")
    shell:
        """
        fasterq-dump {params.outdir} --outdir {params.outdir} \
            --mem {resources.mem_mb}MB --threads {threads} \
            --temp {resources.tmpdir} \
            --force &> {log}; \
        touch {output.flag}
        """

rule fasterq_dump_pe:
    input:
        os.path.join(config['outdir'], "Paired", "{run}", "{run}.sra"),
    output:
        flag = os.path.join(config['outdir'], "Paired", "{run}", "{run}.dumped"),
        fq1 = os.path.join(config['outdir'], "Paired", "{run}", "{run}_1.fastq"),
        fq2 = os.path.join(config['outdir'], "Paired", "{run}", "{run}_2.fastq"),
    params:
        outdir=os.path.join(config['outdir'], "Paired", "{run}"),
        fastqdump = config['fasterq-dump']
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
    threads: 4
    conda:
        "/scicore/home/zavolan/franch0002/myscripts/sra_tools.yaml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "fasterq_dump.log")
    shell:
        """
        mkdir -p {params.outdir}; {params.fastqdump} {params.outdir} --outdir {params.outdir} \
            --mem {resources.mem_mb}MB --threads {threads} \
            --temp {resources.tmpdir} \
            --force &> {log}; \
        touch {output.flag}
        """

#### PIGZ ####
rule compress_fastq_se:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        files = os.path.join(config['outdir'], "Single", "{run}", "{run}.fastq"),
        tmpf = os.path.join(config['outdir'], "Single", "{run}", "{run}.dumped"),
    output:
        flag = os.path.join(config['outdir'], "Single", "{run}", "{run}.pigzflag"),
    threads: 8
    conda:  
        "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "compress_fastq.log")
    shell:
        """
        pigz --best --processes {threads} {input.files}; touch {output.flag} &> {log}
        """

rule compress_fastq_pe:
    "Compress fastq inplace with pigz at best (9) compression level."
    input:
        fq1 = lambda wildcards: getFastq1FromRunList(PE_nested_fastq, wildcards.run),
        fq2 = lambda wildcards: getFastq2FromRunList(PE_nested_fastq, wildcards.run),        
    output:
        flag = os.path.join(config['outdir'], "Paired", "{run}", "{run}.pigzflag"),
    threads: 8
    conda:  
        "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "compress_fastq.log")
    shell:  
        """
        pigz --force --best --processes {threads} {input.fq1}; pigz --force --best --processes {threads} {input.fq2}; touch {output.flag} &> {log}
        """

#### SEQTK ####
rule seqtk_pe:
    input:
        fastq = lambda wildcards: getFromDict(PE_nested_gz, [wildcards.run, wildcards.run_end]),
        tmpf = os.path.join(config['outdir'], "Paired", "{run}", "{run}.pigzflag"),
    output:
        qc = os.path.join(config['outdir'], "Paired", "{run}", "seqtk", "{run_end}.qc")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 1
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "{run_end}.seqtk.log")
    params:
        outdir =  os.path.join(config['outdir'], "Paired", "{run}", "seqtk"),
        seqtk = config['seqtk']
    shell:
        """(mkdir -p {params.outdir}; {params.seqtk} fqchk {input.fastq} > {output.qc} &> {log})"""

rule seqtk_se:
    input:
        fastq = lambda wildcards: SE_fastq.loc[SE_fastq["Run"] == wildcards.run, "Path"],
        tmpf = os.path.join(config['outdir'], "Single", "{run}", "{run}.pigzflag"),
    output:
        qc = os.path.join(config['outdir'], "Single", "{run}", "seqtk", "{run}.qc")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 1
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "seqtk.log"),
    params:
        outdir =  os.path.join(config['outdir'], "Single", "{run}", "seqtk"),
        seqtk = config['seqtk']
    shell:
        """(mkdir -p {params.outdir}; {params.seqtk} fqchk {input.fastq} > {output.qc} &> {log})"""


#### FASTQC ####
rule fastqc_se:
    input:
        fastq = lambda wildcards: SE_fastq.loc[SE_fastq["Run"] == wildcards.run, "Path"]
    output:
        outdir = directory(os.path.join(config['outdir'], "Single", "{run}", "fastqc"))
    params: 
        fastqc = config['fastqc']
    threads: 1
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "fastqc_mapped.log")
    shell:
        """(mkdir -p {output.outdir}; \
        {params.fastqc} --outdir {output.outdir} \
        --threads {threads} \
        {input.fastq}) \
        &> {log}"""

rule fastqc_pe:
    input:
        fastq = lambda wildcards: getFromDict(PE_nested_gz, [wildcards.run, wildcards.run_end])
    output:
        outdir = directory(os.path.join(config['outdir'], "Paired", "{run}", "fastqc", "{run_end}"))
    params: 
        fastqc = config['fastqc']
    threads: 1
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "{run_end}.fastqc_mapped.log")
    shell:
        """(mkdir -p {output.outdir}; \
        {params.fastqc} --outdir {output.outdir} \
        --threads {threads} \
        {input.fastq}) \
        &> {log}"""


#### STAR ####
rule star_gtf_pe:
    input:
        fq1 = lambda wildcards: getFastq1FromRunList(PE_nested_gz, wildcards.run),
        fq2 = lambda wildcards: getFastq2FromRunList(PE_nested_gz, wildcards.run),
        genome_annotated = config['starIndexGTF']
    output:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    params:
        outFileNamePrefix = lambda wildcards, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
        outdir =  os.path.join(config['outdir'],"Paired", "{run}", "BAMS", "withGTF"),
        star = config['star']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "{run}.map_genome_star.log")
    shell:
        """(mkdir -p {params.outdir}; \
        {params.star} \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN {threads} \
        --runMode alignReads \
        --alignEndsType Local \
        --outSAMattrIHstart 0 \
        --outSAMattributes All \
        --runThreadN {threads} \
        --genomeDir {input.genome_annotated} \
        --readFilesIn {input.fq1} {input.fq2} \
        --outFileNamePrefix {params.outFileNamePrefix} &> {log})"""


rule star_nogtf_pe:
    input:
        fq1 = lambda wildcards: getFastq1FromRunList(PE_nested_gz, wildcards.run),
        fq2 = lambda wildcards: getFastq2FromRunList(PE_nested_gz, wildcards.run),
        genome_NOT_annotated = config['starIndexNOGTF']
    output:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    params:
        outFileNamePrefix = lambda wildcards, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
        outdir =  os.path.join(config['outdir'],"Paired", "{run}", "BAMS", "noGTF"),
        star = config['star']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.map_genome_star.log")
    shell:
        """(mkdir -p {params.outdir}; \
        {params.star} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN {threads} \
        --runMode alignReads \
        --alignEndsType Local \
        --outSAMattrIHstart 0 \
        --outSAMattributes All \
        --runThreadN {threads} \
        --genomeDir {input.genome_NOT_annotated} \
        --readFilesIn {input.fq1} {input.fq2} \
        --outFileNamePrefix {params.outFileNamePrefix} &> {log})"""

rule star_gtf_se:
    input:
        fq = lambda wildcards: SE_fastq.loc[SE_fastq["Run"] == wildcards.run, "Path"],
        genome_annotated = config['starIndexGTF']
    output:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    params:
        outFileNamePrefix = lambda wildcards, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
        outdir =  os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF"),
        star = config['star']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "{run}.map_genome_star.log")
    shell:
        """(mkdir -p {params.outdir}; \
        {params.star} \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN {threads} \
        --runMode alignReads \
        --alignEndsType Local \
        --outSAMattrIHstart 0 \
        --outSAMattributes All \
        --runThreadN {threads} \
        --genomeDir {input.genome_annotated} \
        --readFilesIn {input.fq} \
        --outFileNamePrefix {params.outFileNamePrefix} &> {log})"""

rule star_nogtf_se:
    input:
        fq = lambda wildcards: SE_fastq.loc[SE_fastq["Run"] == wildcards.run, "Path"],
        genome_NOT_annotated = config['starIndexNOGTF']
    output:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    params:
        outFileNamePrefix = lambda wildcards, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
        outdir =  os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF"),
        star = config['star']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "{run}.map_genome_star.log")
    shell:
        """(mkdir -p {params.outdir}; \
        {params.star} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN {threads} \
        --runMode alignReads \
        --alignEndsType Local \
        --outSAMattrIHstart 0 \
        --outSAMattributes All \
        --runThreadN {threads} \
        --genomeDir {input.genome_NOT_annotated} \
        --readFilesIn {input.fq} \
        --outFileNamePrefix {params.outFileNamePrefix} &> {log})"""


#### SAMTOOLS ####
## SE ##
rule samtools_index_nogtf_se:
    input:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(config['outdir'], "Single","{run}", "BAMS", "noGTF", "samtools", "{run}.Aligned.sortedByCoord.out.bai")
    params:
        outdir = os.path.join(config['outdir'], "Single","{run}", "BAMS", "noGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_index_nogtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} index {input.bam} -b -o {output.bai} -@ {threads} &> {log})"""

rule samtools_idxstats_nogtf_se:
    input:
        bam = os.path.join(config['outdir'], "Single","{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        stats = os.path.join(config['outdir'], "Single","{run}", "BAMS", "noGTF", "samtools", "{run}.Aligned.sortedByCoord.out.idxstats")
    params:
        outdir = os.path.join(config['outdir'], "Single","{run}", "BAMS", "noGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:1
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_stats_nogtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} idxstats {input.bam} > {output.stats} 2> {log})"""


rule samtools_index_gtf_se:
    input:
        bam = os.path.join(config['outdir'], "Single","{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(config['outdir'], "Single","{run}", "BAMS", "withGTF", "samtools", "{run}.Aligned.sortedByCoord.out.bai")
    params:
        outdir = os.path.join(config['outdir'],"Single", "{run}", "BAMS", "withGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_index_gtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} index {input.bam} -b -o {output.bai} -@ {threads} &> {log})"""

rule samtools_idxstats_gtf_se:
    input:
        bam = os.path.join(config['outdir'], "Single","{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        stats = os.path.join(config['outdir'], "Single","{run}", "BAMS", "withGTF", "samtools", "{run}.Aligned.sortedByCoord.out.idxstats")
    params:
        outdir = os.path.join(config['outdir'], "Single","{run}", "BAMS", "withGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:1
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_stats_gtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} idxstats {input.bam} > {output.stats} 2> {log})"""

## PE ##
rule samtools_index_nogtf_pe:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "samtools", "{run}.Aligned.sortedByCoord.out.bai")
    params:
        outdir = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_index_nogtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} index {input.bam} -b -o {output.bai} -@ {threads} &> {log})"""

rule samtools_idxstats_nogtf_pe:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        stats = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "samtools", "{run}.Aligned.sortedByCoord.out.idxstats")
    params:
        outdir = os.path.join(config['outdir'],"Paired", "{run}", "BAMS", "noGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:1
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_stats_nogtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} idxstats {input.bam} > {output.stats} 2> {log})"""


rule samtools_index_gtf_pe:
    input:
        bam = os.path.join(config['outdir'],"Paired", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        bai = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "samtools", "{run}.Aligned.sortedByCoord.out.bai")
    params:
        outdir = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads: 12
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_index_gtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} index {input.bam} -b -o {output.bai} -@ {threads} &> {log})"""

rule samtools_idxstats_gtf_pe:
    input:
        bam = os.path.join(config['outdir'],"Paired",  "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam")
    output:
        stats = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "samtools", "{run}.Aligned.sortedByCoord.out.idxstats")
    params:
        outdir = os.path.join(config['outdir'],"Paired", "{run}", "BAMS", "withGTF", "samtools"),
        samtools = config['samtools']
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:1
    log: 
        os.path.join(config["log_dir"], "samples", "{run}", "samtools_stats_gtf.log")
    shell:
        """(mkdir -p {params.outdir}; {params.samtools} idxstats {input.bam} > {output.stats} 2> {log})"""


#### MEGADEPTH #### 

## SE ##
# genes #
rule megadepth_gtf_se_genes:
    input:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.genes.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Single",  "{run}", "BAMS", "withGTF", "megadepth", "genes"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "megadepth", "genes", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_genes_gtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

rule megadepth_nogtf_se_genes:
    input:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.genes.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Single",  "{run}", "BAMS", "noGTF", "megadepth", "genes"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "megadepth", "genes", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_genes_nogtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """
# exons #
rule megadepth_gtf_se_exons:
    input:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.exons.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Single",  "{run}", "BAMS", "withGTF", "megadepth", "exons"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "withGTF", "megadepth", "exons", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_exons_gtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

rule megadepth_nogtf_se_exons:
    input:
        bam = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.exons.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Single",  "{run}", "BAMS", "noGTF", "megadepth", "exons"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Single", "{run}", "BAMS", "noGTF", "megadepth", "exons", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_exons_nogtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """
## PE ##
# genes # 
rule megadepth_gtf_pe_genes:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.genes.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF","megadepth", "genes"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "megadepth", "genes", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_genes_gtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

rule megadepth_nogtf_pe_genes:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.genes.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "megadepth", "genes"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "megadepth", "genes",  "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_genes_nogtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

# exons #
rule megadepth_gtf_pe_exons:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.exons.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Paired",  "{run}", "BAMS", "withGTF", "megadepth", "exons"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "withGTF", "megadepth", "exons", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_exons_gtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

rule megadepth_nogtf_pe_exons:
    input:
        bam = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "{run}.Aligned.sortedByCoord.out.bam"),
        bed = config['gtf.exons.bed']
    params:
        prefix = lambda wildcards, output: output.auc.replace("_auc.out", ""),
        outdir = os.path.join(config['outdir'], "Paired",  "{run}", "BAMS", "noGTF", "megadepth", "exons"),
        megadepth = config['megadepth']
    output:
        auc = os.path.join(config['outdir'], "Paired", "{run}", "BAMS", "noGTF", "megadepth", "exons", "{run}_auc.out")
    conda: "/scicore/home/zavolan/franch0002/myscripts/env.yml"
    threads:12
    log:
        os.path.join(config["log_dir"], "samples", "{run}", "megadepth_exons_nogtf.log")
    shell:
        """
        (mkdir -p {params.outdir};
        {params.megadepth} {input.bam} --annotation {input.bed} --threads {threads} --prefix {params.prefix} --auc --frag-dist --include-sofclip --junctions --all-junctions --bigwig --no-index > {output.auc} 2> {log})
        """

## TARGET GENE AND DGE ANALYSIS (on megadepth outputs) ##
rule DGEA:
    input: 
        samples = config['samples'],
        genes_bed = config['gtf.genes.bed'],
    params:
        sample_dir = config['outdir'],
        dea_outdir = os.path.join(config["outdir"], "DEA",),
        script = "/scicore/home/zavolan/franch0002/myscripts/DGEA.R"
    log:
        os.path.join(config['outdir'], "DEA", "DGEA.log"),
    output:
        flag = os.path.join(config['outdir'], "DEA", "prova.flag")
    conda: "/scicore/home/zavolan/franch0002/myscripts/R4pipeline.yaml",
    shell:
        """
        (set +o pipefail;
        mkdir -p {params.dea_outdir};
        Rscript --verbose {params.script} {params.dea_outdir} {input.samples} {input.genes_bed} {params.sample_dir}) &> {log};
        touch {output.flag}
        """ 

