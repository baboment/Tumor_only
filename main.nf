nextflow.enable.dsl=2

// No absolute defaults in public repo
params.samplesheet = null
params.outdir      = 'results'

// References from env OR left null to be set in config/local
params.ref       = System.getenv('REF_FA')        ?: null
params.mills     = System.getenv('MILLS_VCF')     ?: null
params.dbsnp     = System.getenv('DBSNP_VCF')     ?: null
params.gnomad_af = System.getenv('GNOMAD_AF_VCF') ?: null
params.exome_bed = System.getenv('EXOME_BED')     ?: null   // optional
params.pon_vcf   = System.getenv('PON_VCF')       ?: null   // optional

// Parabricks / alignment options
params.num_gpus    = 2
params.bwa_opts    = '-Y'
params.rg_lb       = 'lib1'
params.rg_pl       = 'ILLUMINA'



process FQ2BAM {
  label 'gpu'
  container params.pb_container
  publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram*'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    tuple val(sample), path("${sample}.markdup.cram"), path("${sample}.markdup.cram.crai"), path("${sample}.recal.txt")

  script:
  def intArg = params.exome_bed ? "--interval-file ${params.exome_bed}" : ''
  """
  pbrun fq2bam \
    --num-gpus ${params.num_gpus} \
    --ref ${params.ref} \
    --in-fq ${r1} ${r2} \
    --read-group-id ${sample}.rg1 \
    --read-group-sm ${sample} \
    --read-group-lb ${params.rg_lb} \
    --read-group-pl ${params.rg_pl} \
    --bwa-options "-Y -K 10000000" \
    --out-bam ${sample}.markdup.cram \
    --out-recal-file ${sample}.recal.txt \
    --knownSites ${params.mills} \
    --knownSites ${params.dbsnp} \
    --gpuwrite \
    --gpusort
  """
}

process MUTECTCALLER {
  label 'gpu'
  container params.pb_container
  publishDir "${params.outdir}/vcf_raw", mode: 'copy', pattern: '*.unfiltered.vcf.gz*'

  input:
    val pon_vcf
    tuple val(sample), path(cram), path(crai), path(recal)

  output:
    tuple val(sample), path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), path("${sample}.f1r2.tar.gz"), path("${sample}.unfiltered.vcf.gz.stats")

  script:
  def intArg = params.exome_bed ? "--interval-file ${params.exome_bed}" : ''
  def ponArg = (pon_vcf && pon_vcf != '') ? "--pon ${pon_vcf}" : ''
  """
  pbrun mutectcaller \
    --num-gpus ${params.num_gpus} \
    --ref ${params.ref} \
    --tumor-name ${sample} \
    --in-tumor-bam ${cram} \
    --in-tumor-recal-file ${recal} \
    --mutect-germline-resource ${params.gnomad_af} \
    --out-vcf ${sample}.unfiltered.vcf.gz \
    --mutect-f1r2-tar-gz ${sample}.f1r2.tar.gz \
    ${intArg} \
    ${ponArg}
  """
}

process GET_PILEUPS {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_filtered", mode: 'copy', pattern: '*.pileups.table*'
  input:
    tuple val(sample), path(cram), path(crai), path(recal)

  output:
    tuple val(sample), path("${sample}.pileups.table")

  script:
  def intArg = params.exome_bed ? "-L ${params.exome_bed}" : ''
  """
  gatk GetPileupSummaries \
    -R ${params.ref} \
    -I ${cram} \
    -V ${params.gnomad_af} \
    ${intArg} \
    -O ${sample}.pileups.table
  """
}

process CALC_CONTAM {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_filtered", mode: 'copy', pattern: '*.contamination.table*'
  input:
    tuple val(sample), path(pile)

  output:
    tuple val(sample), path("${sample}.contamination.table")

  script:
  """
  gatk CalculateContamination \
    -I ${pile} \
    -O ${sample}.contamination.table
  """
}

process LEARN_ORIENT {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_filtered", mode: 'copy', pattern: '*.artifacts.tar.gz*'
  input:
    tuple val(sample), path(vcf_in), path(vcf_tbi), path(f1r2), path(stats)

  output:
    tuple val(sample), path("${sample}.artifacts.tar.gz"), path(vcf_in), path(vcf_tbi), path(stats)

  script:
  """
  gatk LearnReadOrientationModel \
    -I ${f1r2} \
    -O ${sample}.artifacts.tar.gz
  """
}

process FILTER_M2 {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_filtered", mode: 'copy', pattern: '*.filtered.vcf.gz*'

  input:
    tuple val(sample), path(artifacts), path(vcf_in), path(vcf_tbi), path(stats), path(contam)

  output:
    tuple val(sample), path("${sample}.filtered.vcf.gz"), path("${sample}.filtered.vcf.gz.tbi")

  script:
  """
  gatk FilterMutectCalls \
    -R ${params.ref} \
    -V ${vcf_in} \
    --contamination-table ${contam} \
    --ob-priors ${artifacts} \
    -O ${sample}.filtered.vcf.gz
  """
}

process SELECT_PASS {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_PASS", mode: 'copy', pattern: '*.pass.vcf.gz*'

  input:
    tuple val(sample), path(filtered), path(filtered_tbi)

  output:
    // SINGLE channel: emit one tuple per sample containing sample + both files
    tuple val(sample), path("${sample}.pass.vcf.gz"), path("${sample}.pass.vcf.gz.tbi")

  script:
  """
  gatk SelectVariants \
    -V ${filtered} \
    --exclude-filtered \
    -O ${sample}.pass.vcf.gz
  """
}

process COPY_REF_TO_CRAM_DIR {
  label 'cpu'
  container params.gatk_container
  publishDir "${params.outdir}/cram", mode: 'copy', overwrite: true

  input:
    // remove `from ref_ch` / `from fq2_done_ch`
    path ref_fa
    val  _

  output:
    path "${ref_fa.getFileName()}"

  script:
  """
  set -euo pipefail
  base=\$(basename "${ref_fa}")
  tmp=".nf_tmp_\$base"

  # Make a *new* file first, then replace the staged input path
  ln -f "${ref_fa}" "\$tmp" 2>/dev/null || cp -a --reflink=auto "${ref_fa}" "\$tmp"

  """
}

// ---------------------- TOP-LEVEL ----------------------
workflow {
  if( !file(params.samplesheet).exists() )
    exit 1, "Samplesheet not found: ${params.samplesheet}"

  reads_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true, sep:',')
    .map { row ->
      if( !row.sample || !row.r1 || !row.r2 )
        throw new IllegalArgumentException('Samplesheet header must be: sample,r1,r2 and cells cannot be empty')
      tuple( row.sample as String, file(row.r1), file(row.r2) )
    }

  pon_vcf_ch = Channel.value( params.pon_vcf ?: '' )
  ref_ch = Channel.of(params.ref)
  fq2   = FQ2BAM(reads_ch)
  fq2_done_ch = fq2.collect().map { true }
  COPY_REF_TO_CRAM_DIR(ref_ch, fq2_done_ch)
  m2    = MUTECTCALLER(pon_vcf_ch, fq2)
  vcf4f = m2
  pile   = GET_PILEUPS(fq2)
  contam = CALC_CONTAM(pile)
  orient = LEARN_ORIENT(vcf4f)

  joined = orient
    .join(contam, by: 0)
    .map { row ->
      // row order after join is flattened: sample, artifacts, vcf_in, vcf_tbi, stats, contamination
      def (s, artifacts, vcf_in, vcf_tbi, stats, cont) = row
      tuple(s, artifacts, vcf_in, vcf_tbi, stats, cont)
    }

  filt  = FILTER_M2(joined)
  passv = SELECT_PASS(filt)

  emit:
    passv
}



