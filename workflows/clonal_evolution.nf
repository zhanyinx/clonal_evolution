

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.date = new java.util.Date().format('yyMMdd')

include {ascat_pyclone; isfile; isdirectory} from '../modules/ascat_pyclone.nf'





workflow CLONAL_EVOLUTION{
    // clonal tmb
    if (isfile(params.input)){
        samples=Channel.from(file(params.input))
    }
    if (isdirectory(params.input)){
        samples=Channel.fromPath("${params.input}/*csv")
    }
    if (params.biomarkers_ascat_config && !params.biomarkers_ascat_config.isEmpty()){
        config=Channel.from(file(params.biomarkers_ascat_config))
    }else{
        config=Channel.from(file("$projectDir/resources/configs/ascat.sarek.config"))
    }
    ascat_pyclone(samples, config)
}