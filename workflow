# create conda environment
conda create -n seq_align -c bioconda
# packages for processing sequencing
# BWA 0.7.17-r1188
# Picard 3.0.0
# star 2.7.10a
# trim-galore 0.6.7
# samtools version 1.6
# subread version 2.0.2
# reference genome (human): http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
# GTF file: http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/ 

##############################################################
#### simple example of workflow for germline variants detection
##############################################################
# using only chromosome 22 (Homo_sapiens.GRCh38.dna.chromosome.22.fa) as example
# "config.yaml" includes sample information

configfile: "config.yaml"

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map:
    input:
        "/home/variants_calling/ref_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.sam"
    shell:
        "bwa mem {input} > {output}"

rule order_coordinate:
    input:
        "mapped_reads/{sample}.sam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "java -jar $PICARD_ROOT/picard.jar SortSam "
        "SORT_ORDER=coordinate "
        "INPUT={input} "
        "OUTPUT={output}"

rule mark_duplications:
    input:
       "sorted_reads/{sample}.bam"
    output:
       "markduplications_reads/{sample}.bam"
    shell:
       "java -jar $PICARD_ROOT/picard.jar MarkDuplicates  "
        "INPUT={input} "
        "OUTPUT={output} "
        "M=marked_dup_metrics.txt "

rule variants_calling:
    input:
       fa="/home/variants_calling/ref_genome/Homo_sapiens.GRCh38.dna.chromosome.22.fa",
       bam=expand("markduplications_reads/{sample}.bam",sample=config["samples"])
    output:
       "calls/all.vcf"
    shell:
       "freebayes -f {input.fa} {input.bam} >{output}"

snakemake --core 1 calls/all.vcf

#######################################
#### simple example of workflow for bulk RNA-seq 
#######################################
# mapping reads for paired-end samples
# config.yaml includes sample information

configfile: "config.yaml"
rule all:
    input: 
        "QC/metrics_summary.py"
# trimming adaptors
rule trim: 
    input:
        expand("{sample}.R{read_no}.fastq.gz", sample=config["samples"], read_no=['1', '2'])
    output:
        trim1out="/trim/{sample}_1_val_1.fq.gz"
        trim2out="/trim/{sample}_2_val_2.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
         
    shell:
        "trim_galore --quality 30"
        "--fastqc"
        "--illumina"
        "--gzip"
        "--length 30"
        "--stringency 3"
        "--retain_unpaired"
        "-o {output}"
        "--paired {input.reads1} {input.reads2}"
        "{log}"

# mapping to reference genom
rule star:
     input:
        fq1="/trim/sample{sample}_1_val_1.fq.gz",
        fq2="/trim/sample{sample}_2_val_2.fq.gz"
        # STAR reference genom index
        idx="index"
        gtf="anno/Homo_sapiens.GRCh38.107.gtf"
     output:
        aln="star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        log="logs/{sample}/Log.out"
        log_final="logs/{sample}/Log.final.out"
     threads: 5
     params:
         extra="--outFileNamePrefix star/{sample}/{sample}_"
     shell:
         "STAR {threads} {params.extra}"
         "--genomeDir {input.idx}"
         "--sjdbGTFfile {input.gtf}"
         "--readFilesIn {input.fq1} {input.fq2}"
         "--outSAMstrandField intronMotif"
         "--readFilesCommand zcat"
         "--outSAMtype BAM SortedByCoordinate"
         "--peOverlapNbasesMin 5"

#  Summrize reads to counts table which is for further analysis
rule featurecounts:
     input:
         aln=expand("star/{sample}/{sample}_Aligned.sortedByCoord.out.bam",sample=SAMPLES)
         gtf="anno/Homo_sapiens.GRCh38.107.gtf"
     output:
         "outputCounts/counts_table.txt"
     params:
         extra="-p --countReadPairs -t exon -g gene_id"
     shell: 
         "featureCounts {params.extra}"
         "a {input.gtf}"
         "-o {output}"
         "{input.aln}"
        
rule mapping_metrics:
     Input:
          "outputCounts/counts_table.txt"
     output:
          "QC/allSamples_mapping_metrics.txt"  
     script:
          "QC/metrics_summary.py"     
         
