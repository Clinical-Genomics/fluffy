params.samplesheet=false
params.fastq=false
params.output=false

if (!params.samplesheet || !params.fastq || !params.output){
    println "missing required parameters"
    println "main.nf --samplesheet <samplesheet_csv> --fastq <fastq_folder> --output <output_folder>"
    exit 0
}

Channel
    .fromPath( params.samplesheet )
    .splitCsv(header: true)
    .map{ row-> tuple(row.FCID, row.SampleID, row.index, row.index2, row.SampleName, row.Project, row.Library_nM, row.SequencingDate) }
    .set { samplesheet_ch }

Channel
    .fromPath( params.samplesheet )
    .splitCsv(header: true)
    .map{ row-> row.SampleID }
    .unique()
    .set { sample_ch }
    
process bwa_aln{

    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
       val(sampleID) from sample_ch

    output:
       val(sampleID),file("${sampleID}.bam"), file("${sampleID}.bam.bai") into bwa_aln_ch

    script:
       def R1 = "<( zcat ${params.fastq}/${sampleID}/*${sampleID}*_R1*fastq.gz )"
       def R2 = "<( zcat ${params.fastq}/${sampleID}/*${sampleID}*_R2*fastq.gz )"

       """
       bwa aln -n 0 -k 0 -t ${task.cpus} ${params.ref} ${R1} > ${sampleID}_R1.sai
       bwa aln -n 0 -k 0 -t ${task.cpus} ${params.ref} ${R2} > ${sampleID}_R2.sai
       bwa sampe -n -1 ${params.ref} ${sampleID}_R1.sai ${sampleID}_R2.sai ${R1} ${R2} | samtools sort --output-fmt BAM -@ ${task.cpus} - > ${sampleID}.bam
       samtools index ${sampleID}.bam
       """

}

