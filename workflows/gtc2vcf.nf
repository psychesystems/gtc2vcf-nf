// Call idat genotypes and convert to VCF
// https://github.com/freeseek/gtc2vcf

params.command = "identify"
params.dataset = "Dataset"
params.idats = "data/idats/*.idat"
params.egt = "data/manifest/GSA.egt"
params.bpm = "data/manifest/GSA.bpm"
params.csv = "data/manifest/GSA.csv"

params.gtc = "data/gtc/GSA/*.gtc"
params.fasta = "data/fasta/ref.fasta"

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

        EGT_CH = Channel.fromPath(params.egt, checkIfExists: true)
        BPM_CH = Channel.fromPath(params.bpm, checkIfExists: true)

        IDAT_CH = Channel.fromFilePairs(params.idats) { file -> file.simpleName.split("_")[0,1] }

        IDAT_MANIFEST_CH = IDAT_CH
          .combine(EGT_CH)
          .combine(BPM_CH)

        GTC_CH = GTC(IDAT_MANIFEST_CH)
        
    }

    if(params.command == "vcf") {

      DATASET_CH = Channel.of(params.dataset)

      GTC_CH = Channel.fromPath(params.gtc)
        .collect()
        .map { it -> [it] }

      FASTA_GZ_CH = Channel.fromPath(params.fasta, checkIfExists: true)

      FASTA_CH = FASTA(FASTA_GZ_CH)

      EGT_CH = Channel.fromPath(params.egt, checkIfExists: true)
      BPM_CH = Channel.fromPath(params.bpm, checkIfExists: true)
      CSV_CH = Channel.fromPath(params.csv, checkIfExists: true)

      GTC_FASTA_CH = DATASET_CH
        .combine(GTC_CH)
        .combine(EGT_CH)
        .combine(BPM_CH)
        .combine(CSV_CH)
        .combine(FASTA_CH)

      VCF_CH = VCF(GTC_FASTA_CH)

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

// idat2gtc
process GTC {

    publishDir "data/gtc/${bpm.simpleName}"

    cpus = 1
    memory = 1.Gb

    input:
    tuple val(sentrix), path(idats), path(egt), path(bpm)

    output:
    path("*.gtc")

    script:
    """
    bcftools +idat2gtc \
      --bpm ${bpm} \
      --egt ${egt} \
      --idats . \
      --output .
    """
}

// preprocess FASTA
process FASTA {
  tag "${fasta.baseName}"

  cpus = 1
  memory = 1.Gb

  input:
  path(fasta)

  output:
  tuple path(fasta, includeInputs: true), path("*.{amb,ann,fai,pac}")

  script:
  """
  samtools faidx ${fasta}
  bwa index ${fasta}
  """
}

// convert to vcf
process VCF {
  tag "${fasta.baseName}"

  publishDir "data/vcf/${dataset}/${fasta.baseName}"

  cpus = 8
  memory = 24.Gb

  input:
  tuple val(dataset), path(gtc), path(egt), path(bpm), path(csv), path(fasta), path(index)

  output:
  path("*.bcf")

  script:
  """
  bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm ${bpm} \
  --csv ${csv} \
  --egt ${egt} \
  --gtcs . \
  --fasta-ref ${fasta} \
  --extra ${dataset}.tsv \ | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm \
  --no-version \
  --output ${dataset}.bcf \
  --output-type b \
  --check-ref x \
  --fasta-ref ${fasta} \
  --write-index \
  --threads ${task.cpus}
  """
}