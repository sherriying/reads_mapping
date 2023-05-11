# create conda environment
conda create -n seq_align -c bioconda
# packages for processing sequencing
# star 2.7.10a
# trim-galore 0.6.7
# samtools version 1.6
# subread version 2.0.2
# reference genome (human): http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
# GTF file: http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/ 

#######################################
#### simple example of workflow #######
#######################################
# mapping reads for paired-end samples
# config.yaml include sample names

configfile: "config.yaml"
SAMPLES=config["samples"]
rule all:
    input: 
        "QC/metrics_summary.py",
        expand("{sample}.R{read_no}.fastq.gz", sample=SAMPLES, read_no=['1', '2'])
# trimming adaptors
rule trim: 
    input:
        reads1="/reads/{sample}.R1.fq.gz",
        reads2="/reads/{sample}.R2.fq.gz"
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
         
