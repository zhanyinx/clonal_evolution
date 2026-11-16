
process generate_pyclone{
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 2.GB * task.attempt }
    time { 4.h * task.attempt }
    publishDir "${params.outdir}/${patient}/pyclone", mode: "copy"
    fair true

    input:
        tuple val(patient), val(sample_ids), val(sex), val(mafs), val(cellularity), path(crams), path(crais), path(pluriploidy), path(cnvs)
    output:
       tuple val(patient), path("${patient}.pyclone.tsv"), emit: input
       path("${patient}.pyclone.input_report.tsv"), emit: report
    script:
    // groupTuple gives lists; a single staged file may arrive as a bare path
    def ids = [sample_ids].flatten()
    def sexes = [sex].flatten()
    def maf_list = [mafs].flatten()
    def cram_list = [crams].flatten()
    def pp_list = [pluriploidy].flatten()
    def cnv_list = [cnvs].flatten()
    def rows = ['sample_id\tmaf\talignment\tcnvs\tpurityploidy\tsex']
    ids.eachWithIndex { id, i ->
        rows << [id, maf_list[i], cram_list[i].name, cnv_list[i].name, pp_list[i].name, sexes[i]].join('\t')
    }
    """
    cat > samples.tsv <<'END_OF_SAMPLES'
${rows.join('\n')}
END_OF_SAMPLES
    create_input4pyclone.py \
        --samples samples.tsv \
        --fasta ${params.fasta} \
        --output ${patient}.pyclone.tsv \
        --report ${patient}.pyclone.input_report.tsv \
        --max_gap ${params.cn_max_gap}
    """
}


process pyclone{
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    time { 4.h * task.attempt }
    publishDir "${params.outdir}/${patient}/pyclone", mode: "copy"
    input:
        tuple val(patient), path(pyclone)
    output:
       tuple val(patient), path("*.pyclone.output.tsv")
    script:
    """
        # run pyclone
        pyclone-vi fit -i ${pyclone} -o ${patient}.pyclone.h5 -c ${params.num_clusters} -d beta-binomial -r ${params.num_restarts}
        pyclone-vi write-results-file -i ${patient}.pyclone.h5 -o ${patient}.pyclone.output.tsv
    """
}
