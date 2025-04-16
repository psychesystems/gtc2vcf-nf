// Call idat genotypes and convert to VCF
// https://github.com/freeseek/gtc2vcf

params.command = "identify"
params.idats = "data/idats/*.idat"

workflow {
    if(params.command == "identify") {
        IDAT_CH = Channel.fromFilePairs(params.idats) { file -> file.simpleName.split("_")[0,1] }

        ID_CH = IDENTIFY(IDAT_CH)
    }
}

process IDENTIFY {
    tag "${sentrix[0]}_${sentrix[1]}"

    cpus = 1
    memory = 1.Gb

    input:
    tuple val(sentrix), path(idats)

    output:
    tuple val(sentrix), path("${sentrix.join('_')}.tsv")

    script:
    """
    bcftools +gtc2vcf \
    -i -g ./ > ${sentrix.join('_')}.tsv
    """
}