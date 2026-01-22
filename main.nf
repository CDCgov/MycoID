#!/usr/bin/env nextflow

params.schema_path = "${workflow.projectDir}/nextflow_schema.json"


process concatenateFastq {

    tag { sample }
    errorStrategy 'ignore'

    publishDir "${params.output}/concatenated", mode: 'copy', pattern: '*.gz'

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("${sample}_combined.fastq.gz")
 
    script:
    """
    cat ${files} > ${sample}_combined.fastq.gz
    """
}
 
process fastp {

    tag { sample }
    errorStrategy 'ignore'

    publishDir "${params.output}/cleaned", mode: 'copy', pattern: '*.gz'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_long_clean.fastq.gz"), emit: out
    tuple val(sample), path("${sample}.json"), emit: json

    script:
    """
    fastplong -i ${fastq} -q 30 --length_required ${params.length_required} --length_limit ${params.length_limit} -o ${sample}_long_clean.fastq.gz -j ${sample}.json
    """
}

process downsample {

    tag { sample }
    errorStrategy 'ignore'

    publishDir "${params.output}/downsampled", mode: 'copy', pattern: '*.fastq'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_downsampled.fastq")

    script:
    """
    ontime --to 12h -o ${sample}_downsampled.fastq ${fastq}
    """
}

process consensus {

    tag { sample }
    errorStrategy 'ignore'

    publishDir "${params.output}/consensus", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.fasta")

    script:
    """ 
    NGSpeciesID --ont --consensus --medaka --fastq ${fastq} --outfolder ${sample}
    cat ${sample}/*.fasta > ${sample}.fasta
    """

}

process blast {

    tag { sample }

    maxRetries 3
    errorStrategy 'ignore'

    publishDir "${params.output}/blast", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_classification.csv")

    script:
    """
    # run blast
    update_blastdb.pl --decompress taxdb
    blastn -query ${fasta} -db core_nt -entrez_query "Fungi[Organism]" -remote -evalue 0.00001 -outfmt "10 sscinames sseqid staxids evalue qseq length pident slen" > ${sample}_blast.csv
    
    # filter, sort and format the output
    awk -F, '\$1 !~ /uncultured|sp\\.|fungal|fungus|subsp\\./ && \$7 >= ${params.percent} && \$6 >= 0.8*\$8' ${sample}_blast.csv | \
    sort -t',' -k7,7nr -k4,4n -k6,6nr | \
    cut -d',' -f1,2,4-7 | \
    awk -v sample="${sample}" '{print sample "," \$0}' > ${sample}_classification.csv
    """

}

process sample_report {

    tag { sample }
    errorStrategy 'ignore'

    publishDir "${params.output}/report/sample", mode: 'copy'

    input:
    tuple val(sample), path(blast)

    output:
    path("${sample}_summary.csv")

    script:
    def scriptName = "MycoID - Fungal ID Analysis"
    def user = params.user
    def version = params.version
    def runDate = new Date().format('yyyy-MM-dd')
    """
    touch ${sample}_summary.csv
    echo "sscinames,sseqid,evalue,qseq,length,pident" > ${sample}_summary.csv
    cat ${blast} | cut -d',' -f2- >> ${sample}_summary.csv
    echo -e "${scriptName}\nUser: ${user}\nVersion: ${version}\nDate: ${runDate}\nSample: ${sample}\n"  | cat - ${sample}_summary.csv > temp.txt && mv temp.txt ${sample}_summary.csv
    """
}


process combined_report {

    errorStrategy 'ignore'

    publishDir "${params.output}/report/combined", mode: 'copy'

    input:
    path(blast)

    output:
    path("final_summary.csv")

    script:
    def scriptName = "MycoID - Fungal ID Analysis"
    def user = params.user
    def version = params.version
    def runDate = new Date().format('yyyy-MM-dd')
    """
    #combined
    touch final_summary.csv
    echo "sscinames,sseqid,evalue,qseq,length,pident" > final_summary.csv
    cat ${blast} >> final_summary.csv
    echo -e "${scriptName}\nUser: ${user}\nVersion: ${version}\nDate: ${runDate}\n"  | cat - final_summary.csv > temp.txt && mv temp.txt final_summary.csv
    """

}

process qcReport {

    tag { sample }

    publishDir "${params.output}/qc_report", mode: 'copy'

    input:
    tuple val(sample), path(json)

    output:
    path("${sample}_qc.txt")

    script:
    """
    jq -r '
    "summary\\n",
    "before_filtering:",
    "total_reads:\\(.summary.before_filtering.total_reads)",
    "total_bases:\\(.summary.before_filtering.total_bases)",
    "q20_bases:\\(.summary.before_filtering.q20_bases)",
    "q30_bases:\\(.summary.before_filtering.q30_bases)",
    "q20_rate:\\(.summary.before_filtering.q20_rate | tostring | .[0:7])",
    "q30_rate:\\(.summary.before_filtering.q30_rate | tostring | .[0:7])",
    "read_mean_length:\\(.summary.before_filtering.read_mean_length)",
    "gc_content:\\(.summary.before_filtering.gc_content | tostring | .[0:7])",
    ' ${json} > ${sample}_qc.txt
    """

}


workflow {

    if (!params.input) {
        error "ERROR: Missing required input parameter. Please specify the input directory using '--input'."
    }
    if (!params.output) {
        error "ERROR: Missing required output parameter. Please specify the output directory using '--output'."
    }
    if (!params.user) {
        error "ERROR: Missing required 'User' argument. Please specify your CDC USER ID using '--user'."
    }

    grouped_samples = Channel.fromPath("${params.input}/*/*.fastq.gz", checkIfExists:true) \
    | map { file -> 
      def key = file.parent.name
      return tuple(key, file)
    } \
    | groupTuple() 

    concatenated = concatenateFastq(grouped_samples)
    cleaned = fastp(concatenated)
    downsampled = downsample(cleaned.out)
    assemblies = consensus(downsampled)
    blastOut = blast(assemblies)

    // Reports
    sample_report(blastOut)
    combined_report(blastOut.map { it[1] }.collect())
}
