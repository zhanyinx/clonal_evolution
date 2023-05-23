

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.date = new java.util.Date().format('yyMMdd')

include {ascat_pyclone;} from '../modules/ascat_pyclone.nf'


workflow CLONAL_EVOLUTION{
    // clonal tmb
    sample=Channel.from(file(params.input))
    if (params.biomarkers_ascat_config && !params.biomarkers_ascat_config.isEmpty()){
        config=Channel.from(file(params.biomarkers_ascat_config))
    }else{
        config=Channel.from(file("$projectDir/resources/configs/ascat.sarek.config"))
    }
    ascat_pyclone(sample, config)
}