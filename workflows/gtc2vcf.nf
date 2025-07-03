// Call idat genotypes and convert to VCF
// https://github.com/freeseek/gtc2vcf

params.command = "identify"
params.idats = "data/idats/*.idat"
params.egt = "https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/InfiniumGlobalScreeningArrayv4.0ClusterFile-egt.zip"
params.bpm = "https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium-gsa-with-gcra/GSA-48v4-0_20085471_D2.bpm"


workflow {
    if(params.command == "identify") {
        IDAT_CH = Channel.fromFilePairs(params.idats) { file -> file.simpleName.split("_")[0,1] }

        ID_CH = IDENTIFY(IDAT_CH)
        
        ID_CH
          .splitCsv( header: true, sep: "\t")
          .map {it -> [it[0], it[1].chip_type_guess] }
          .groupTuple()
          .map { it -> [it[1], it[0]] }
          .groupTuple()
          .map { it -> "Chip type guess (G/R): ${it[0]} Count: ${it[1].size()}"}
          .view()
    }

    if(params.command == "gtc") {
        IDAT_CH = Channel.fromFilePairs(params.idats) { file -> file.simpleName.split("_")[0,1] }

        // manifest files
        EGT_CH = EGT(Channel.fromPath(params.egt))
        BPM_CH = Channel.fromPath(params.bpm)

        IDAT_MANIFEST_CH = IDAT_CH
          .combine(EGT_CH)
          .combine(BPM_CH)

        GTC_CH = GTC(IDAT_MANIFEST_CH)
        
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

// unzip EGT
process EGT {

    cpus = 1
    memory = 1.Gb

    input:
    path(zip)

    output:
    path("*.egt")

    script:
    """
    unzip ${zip}
    """
}

// idat2gtc
process GTC {

    cpus = 1
    memory = 1.Gb

    input:
    tuple val(sentrix), path(idats), path(egt), path(bpm)

    output:
    tuple path(egt), path(bpm), path("*.gtc")

    script:
    """
    bcftools +idat2gtc \
      --bpm ${bpm} \
      --egt ${egt} \
      --idats . \
      --output .
    """
}