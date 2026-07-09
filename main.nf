#!/usr/bin/env nextflow

params.schema_path = "${workflow.projectDir}/nextflow_schema.json"

// ─── PROCESSES ───────────────────────────────────────────────────────────────

process concatenateFastq {

    tag "${sample}"
    publishDir "${params.output}/concatenated", mode: 'copy', pattern: '*.gz'

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("${sample}_combined.fastq.gz"), env(STATUS), emit: out

    script:
    """
    if cat ${files} > ${sample}_combined.fastq.gz || true; then
        if [ -s ${sample}_combined.fastq.gz ]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
        fi
    fi
    """
}

process fastp {

    tag "${sample}"
    publishDir "${params.output}/cleaned", mode: 'copy', pattern: '*.gz'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_long_clean.fastq.gz"), env(STATUS), emit: out
    tuple val(sample), path("${sample}.json"), emit: json

    script:
    """
    if [ ! -s ${fastq} ]; then
        STATUS="FAIL"
        touch ${sample}_long_clean.fastq.gz
        echo '{}' > ${sample}.json
    elif fastplong -i ${fastq} -q 30 --length_required ${params.length_required} --length_limit ${params.length_limit} -o ${sample}_long_clean.fastq.gz -j ${sample}.json || true; then
        if [ -s ${sample}_long_clean.fastq.gz ]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
            touch ${sample}_long_clean.fastq.gz
            echo '{}' > ${sample}.json
        fi
    fi
    """
}

process downsample {

    tag "${sample}"
    publishDir "${params.output}/downsampled", mode: 'copy', pattern: '*.fastq'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}_downsampled.fastq"), env(STATUS), emit: out

    script:
    """
    if [ ! -s ${fastq} ]; then
        STATUS="FAIL"
        touch ${sample}_downsampled.fastq
    elif ontime --to 12h -o ${sample}_downsampled.fastq ${fastq} || true; then
        if [ -s ${sample}_downsampled.fastq ]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
            touch ${sample}_downsampled.fastq
        fi
    fi
    """
}

process consensus {

    tag "${sample}"
    publishDir "${params.output}/consensus", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.fasta"), env(STATUS), emit: out

    script:
    """
    if [ ! -s ${fastq} ]; then
        STATUS="FAIL"
        touch ${sample}.fasta
    elif NGSpeciesID --ont --consensus --medaka --fastq ${fastq} --outfolder ${sample} || true; then
        cat ${sample}/*.fasta > ${sample}.fasta 2>/dev/null || true
        if [ -s ${sample}.fasta ]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
            touch ${sample}.fasta
        fi
    fi
    """
}

process blast {

    tag "${sample}"
    publishDir "${params.output}/blast", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_classification.csv"), env(STATUS), emit: out

    script:
    """
    update_blastdb.pl --decompress taxdb

    if [ ! -s ${fasta} ]; then
        STATUS="FAIL"
        touch ${sample}_classification.csv
    elif blastn -query ${fasta} -db core_nt -entrez_query "Fungi[Organism]" -remote -evalue 0.00001 \
        -outfmt "10 sscinames sseqid staxids evalue qseq length pident slen" > ${sample}_blast.csv || true; then

        awk -F, '\$1 !~ /uncultured|sp\\.|fungal|fungus|subsp\\./ && \$7 >= ${params.percent} && \$6 >= 0.8*\$8' ${sample}_blast.csv | \
            sort -t',' -k7,7nr -k4,4n -k6,6nr | \
            cut -d',' -f1,2,4-7 | \
            awk -v sample="${sample}" '{print sample "," \$0}' > ${sample}_classification.csv

        if [ -s ${sample}_classification.csv ]; then
            STATUS="PASS"
        else
            STATUS="FAIL"
            touch ${sample}_classification.csv
        fi
    fi
    """
}

process sample_report {

    tag "${sample}"
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
    echo -e "${scriptName}\nUser: ${user}\nVersion: ${version}\nDate: ${runDate}\nSample: ${sample}\n" | cat - ${sample}_summary.csv > temp.txt && mv temp.txt ${sample}_summary.csv
    """
}

process combined_report {

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
    touch final_summary.csv
    echo "sscinames,sseqid,evalue,qseq,length,pident" > final_summary.csv
    cat ${blast} >> final_summary.csv
    echo -e "${scriptName}\nUser: ${user}\nVersion: ${version}\nDate: ${runDate}\n" | cat - final_summary.csv > temp.txt && mv temp.txt final_summary.csv
    """
}

process qcReport {

    tag "${sample}"
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
    "gc_content:\\(.summary.before_filtering.gc_content | tostring | .[0:7])"
    ' ${json} > ${sample}_qc.txt
    """
}

process failureReport {

    publishDir "${params.output}/report", mode: 'copy'

    input:
    val(failures)

    output:
    path("failed_samples.tsv")

    exec:
    def file = task.workDir.resolve("failed_samples.tsv")
    if (failures == ["No failures"]) {
        file.text = "All samples completed successfully\n"
    } else {
        file.text = "sample\tstage\n"
        failures.each { line -> file.append(line + "\n") }
    }
}

// ─── WORKFLOW ─────────────────────────────────────────────────────────────────

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

    grouped_samples = Channel.fromPath("${params.input}/*/*.fastq.gz", checkIfExists: true)
        | map { file ->
            def key = file.parent.name
            return tuple(key, file)
        }
        | groupTuple()

    // ── concatenate ──
    concatenated = concatenateFastq(grouped_samples)

    concatenated.out.branch {
        pass: it[2] == "PASS"
        fail: it[2] == "FAIL"
    }.set { concat_branch }

    concat_failures = concat_branch.fail.map { sample, file, status -> "${sample}\tconcatenate" }

    // ── fastp ──
    cleaned = fastp(concat_branch.pass.map { sample, file, status -> tuple(sample, file) })

    cleaned.out.branch {
        pass: it[2] == "PASS"
        fail: it[2] == "FAIL"
    }.set { fastp_branch }

    fastp_failures = fastp_branch.fail.map { sample, file, status -> "${sample}\tfastp" }

    // ── downsample ──
    downsampled = downsample(fastp_branch.pass.map { sample, file, status -> tuple(sample, file) })

    downsampled.out.branch {
        pass: it[2] == "PASS"
        fail: it[2] == "FAIL"
    }.set { downsample_branch }

    downsample_failures = downsample_branch.fail.map { sample, file, status -> "${sample}\tdownsample" }

    // ── consensus ──
    assemblies = consensus(downsample_branch.pass.map { sample, file, status -> tuple(sample, file) })

    assemblies.out.branch {
        pass: it[2] == "PASS"
        fail: it[2] == "FAIL"
    }.set { consensus_branch }

    consensus_failures = consensus_branch.fail.map { sample, file, status -> "${sample}\tconsensus" }

    // ── blast ──
    blastOut = blast(consensus_branch.pass.map { sample, file, status -> tuple(sample, file) })

    blastOut.out.branch {
        pass: it[2] == "PASS"
        fail: it[2] == "FAIL"
    }.set { blast_branch }

    blast_failures = blast_branch.fail.map { sample, file, status -> "${sample}\tblast" }

    // ── reports ──
    sample_report(blast_branch.pass.map { sample, file, status -> tuple(sample, file) })
    combined_report(blast_branch.pass.map { it[1] }.collect())
    qcReport(cleaned.json)

    // ── failure report ──
    all_failures = concat_failures
        .mix(fastp_failures)
        .mix(downsample_failures)
        .mix(consensus_failures)
        .mix(blast_failures)
        .collect()
        .ifEmpty(["No failures"])

    failureReport(all_failures)
}
