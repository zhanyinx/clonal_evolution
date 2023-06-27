process timescape{
    cpus '1'
    errorStrategy 'retry'
    maxRetries = 3
    memory { 1.GB * task.attempt }
    publishDir "${params.outdir}/${patient}/timescape", mode: "copy"
    input:
        tuple val(patient), path(citup)
        file(mut_sample_clus)
        tuple val(patient), path(pyclone)
    output:
        tuple val(patient), path("*.html")
    script:
    """
        optimal="\$(citup_optimal.py -i $citup)"
        citup2timescape.R -c $citup -m ${mut_sample_clus} -b \$optimal -s $params.timescape_sorting -o ${patient}.timescape -p $pyclone
    """
}