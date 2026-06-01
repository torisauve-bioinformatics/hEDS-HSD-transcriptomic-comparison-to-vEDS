nextflow.enable.dsl=2

params.fastq = "/scratch/vks266/transcrip_final/heds/output/fastq/*.fastq.gz"
params.star_index = "/scratch/vks266/transcrip_final/heds/ampliseq_pipe/hg19_star_index"
params.gtf = "/scratch/vks266/transcrip_final/heds/ampliseq_pipe/hg19.refGene.gtf"
params.outdir = "/scratch/vks266/transcrip_final/heds/ampliseq_pipe/results2"

process TRIM {
    container 'quay.io/biocontainers/cutadapt:5.2--py313h8c92656_1'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    path fastq

    output:
    path "trimmed_${fastq}"

    script:
    """
    cutadapt \
        -a ATCACCGACTGCCCATAGAGAGG \
        --minimum-length 20 \
        -j ${task.cpus} \
        -o trimmed_${fastq} \
        ${fastq}
    """
}

process ALIGN {
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_8'
    publishDir "${params.outdir}/bam", mode: 'copy'
    input:
    path fastq
    path star_index
    output:
    path "*.sortedByCoord.out.bam"

    script:
    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${fastq} \
        --readFilesCommand zcat \
        --alignEndsType Local \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        --alignIntronMax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outFileNamePrefix ${fastq.simpleName}_
    """
}

process COUNT {
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'
    publishDir "${params.outdir}/counts", mode: 'copy'	
    input:
    path bam
    path gtf
    output:
    path "counts_${bam.baseName}.txt*"

    script:
    def sample = bam.name.replace('.sortedByCoord.out.bam', '')
    """
    featureCounts \
        -T ${task.cpus} \
        -s 0 \
        -O \
        -M \
        --fraction \
        -t exon \
        -g gene_id \
        -a ${gtf} \
        -o counts_${bam.baseName}.txt \
        ${bam}
    """
}

workflow {
    fastqs = Channel.fromPath(params.fastq)
    star_index = Channel.fromPath(params.star_index, type: 'dir').first()
    gtf = Channel.fromPath(params.gtf).first() 
    trimmed = TRIM(fastqs)
    bams = ALIGN(trimmed, star_index)
    COUNT(bams, gtf)
}
