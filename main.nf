nextflow.enable.dsl=2

// Site-specific defaults should live in nextflow.config or a local config file.
params.samplesheet = System.getenv('SAMPLESHEET') ?: null
params.outdir      = 'results'

// References come from env, nextflow.config, params files, or local config.
// Keep Mutect resources aligned to Broad hg38/assembly38, not GRCh38.p14.
params.ref                       = System.getenv('REF_FA')                       ?: null
params.mills                     = System.getenv('MILLS_VCF')                    ?: null
params.dbsnp                     = System.getenv('DBSNP_VCF')                    ?: null
params.mutect_germline_resource  = System.getenv('MUTECT_GERMLINE_RESOURCE_VCF') ?: System.getenv('GNOMAD_AF_VCF') ?: null
params.contamination_sites_vcf   = System.getenv('CONTAMINATION_SITES_VCF')      ?: null
params.exome_bed                 = System.getenv('EXOME_BED')                    ?: null
// If no vendor bait BED is available, fall back to the target BED. This makes
// CollectHsMetrics target-depth outputs useful, but bait/enrichment metrics are
// only approximate until a true bait BED is provided.
params.exome_bait_bed            = System.getenv('EXOME_BAIT_BED')               ?: params.exome_bed
// Fast mode for target-coverage QC: keep the Picard hs metrics used by MultiQC,
// but do not re-scan the CRAM for duplicate samtools summaries unless requested.
params.target_cov_emit_aux       = false
// Picard per-target coverage is substantially slower on large target sets.
params.target_cov_per_target     = false
params.mutect_pon_vcf            = System.getenv('MUTECT_PON_VCF')               ?: System.getenv('PON_VCF') ?: null
params.pon_vcf                   = System.getenv('PON_VCF')                      ?: params.mutect_pon_vcf
params.mutect_extra_args         = ''
params.deepsomatic_model_type    = 'WES_TUMOR_ONLY'
params.deepsomatic_fallback_model_type = 'WGS_TUMOR_ONLY'
params.deepsomatic_extra_args    = ''
params.clairs_to_platform        = 'ilmn_ssrs'
params.clairs_to_extra_args      = ''
params.snpeff_db                 = 'hg38'
params.snpeff_jar                = ''
params.snpeff_java_mem           = '8g'
params.snpeff_extra_args         = ''
params.annotation_bcftools_args  = ''
params.python_bin                = 'python3'
params.apptainer_bind_args       = System.getenv('APPTAINER_BIND_ARGS') ?: ''
// Parabricks / alignment options
params.num_gpus    = 2
params.bwa_opts    = '-Y'
params.rg_lb       = 'lib1'
params.rg_pl       = 'ILLUMINA'

def isBlank(value) {
  return value == null || value.toString().trim() == ''
}

def pathExists(value) {
  return !isBlank(value) && file(value.toString()).exists()
}

def presetDefaults = [
  default : [
    mutect_low_lod              : false,
    clairs_to_snv_min_af        : '0.02',
    clairs_to_indel_min_af      : '0.05',
    clairs_to_min_coverage      : '10',
    postfilter_snv_min_af       : '0.02',
    postfilter_indel_min_af     : '0.05',
    postfilter_snv_min_alt_reads: '5',
    postfilter_indel_min_alt_reads: '5',
    postfilter_snv_min_dp       : '30',
    postfilter_indel_min_dp     : '30',
  ],
  moderate: [
    mutect_low_lod              : false,
    clairs_to_snv_min_af        : '0.02',
    clairs_to_indel_min_af      : '0.05',
    clairs_to_min_coverage      : '10',
    postfilter_snv_min_af       : '0.02',
    postfilter_indel_min_af     : '0.05',
    postfilter_snv_min_alt_reads: '5',
    postfilter_indel_min_alt_reads: '5',
    postfilter_snv_min_dp       : '30',
    postfilter_indel_min_dp     : '30',
  ],
  low_af  : [
    mutect_low_lod              : true,
    clairs_to_snv_min_af        : '0.01',
    clairs_to_indel_min_af      : '0.02',
    clairs_to_min_coverage      : '10',
    postfilter_snv_min_af       : '0.01',
    postfilter_indel_min_af     : '0.02',
    postfilter_snv_min_alt_reads: '10',
    postfilter_indel_min_alt_reads: '10',
    postfilter_snv_min_dp       : '50',
    postfilter_indel_min_dp     : '50',
  ],
]

def resolvedPipelinePreset = (params.pipeline_preset ?: 'default').toString().trim()
if (!presetDefaults.containsKey(resolvedPipelinePreset)) {
  exit 1, "Unsupported --pipeline_preset '${resolvedPipelinePreset}'. Choose one of: ${presetDefaults.keySet().join(', ')}"
}

def presetOrParam = { String key ->
  !isBlank(params[key]) ? params[key] : presetDefaults[resolvedPipelinePreset][key]
}

def resolvedMutectLowLod = presetOrParam('mutect_low_lod').toString().toBoolean()
def resolvedClairsToSnvMinAf = presetOrParam('clairs_to_snv_min_af').toString()
def resolvedClairsToIndelMinAf = presetOrParam('clairs_to_indel_min_af').toString()
def resolvedClairsToMinCoverage = presetOrParam('clairs_to_min_coverage').toString()
def resolvedPostfilterSnvMinAf = presetOrParam('postfilter_snv_min_af').toString()
def resolvedPostfilterIndelMinAf = presetOrParam('postfilter_indel_min_af').toString()
def resolvedPostfilterSnvMinAltReads = presetOrParam('postfilter_snv_min_alt_reads').toString()
def resolvedPostfilterIndelMinAltReads = presetOrParam('postfilter_indel_min_alt_reads').toString()
def resolvedPostfilterSnvMinDp = presetOrParam('postfilter_snv_min_dp').toString()
def resolvedPostfilterIndelMinDp = presetOrParam('postfilter_indel_min_dp').toString()
def resolvedMergeMinCallers = (params.merge_min_callers ?: 2) as Integer
if (resolvedMergeMinCallers < 1 || resolvedMergeMinCallers > 3) {
  exit 1, "--merge_min_callers must be between 1 and 3 for the three-caller comparison; got ${resolvedMergeMinCallers}"
}

def clairsCustomPonFiles = [
  params.clairs_to_gnomad,
  params.clairs_to_dbsnp,
  params.clairs_to_pon,
  params.clairs_to_colorsdb,
]
if (clairsCustomPonFiles.any { !isBlank(it) } && !clairsCustomPonFiles.every { pathExists(it) }) {
  log.warn("ClairS-TO custom PoN resources are incomplete or missing; the workflow will omit --panel_of_normals and use ClairS-TO defaults.")
}

def samplesheetValue(row, String name) {
  return row.containsKey(name) && row[name] ? row[name].toString().trim() : null
}

def failSamplesheet(int rowNumber, String message) {
  exit 1, "Samplesheet row ${rowNumber}: ${message}"
}

def requireSamplesheetPath(int rowNumber, String sample, String label, String pathText) {
  if (!pathText) {
    failSamplesheet(rowNumber, "missing ${label} for sample '${sample}'")
  }
  if (!file(pathText).exists()) {
    failSamplesheet(rowNumber, "path for ${label} does not exist for sample '${sample}': ${pathText}")
  }
}

def requireVcfIndexIfCompressed(int rowNumber, String sample, String label, String pathText) {
  requireSamplesheetPath(rowNumber, sample, label, pathText)
  if (pathText.endsWith('.vcf.gz') || pathText.endsWith('.bcf')) {
    def hasIndex = file("${pathText}.tbi").exists() || file("${pathText}.csi").exists()
    if (!hasIndex) {
      failSamplesheet(rowNumber, "missing index for ${label} for sample '${sample}'; expected ${pathText}.tbi or ${pathText}.csi")
    }
  }
}

def parseSamplesheetRow(row, int rowNumber, Set seenSamples) {
  def sample = samplesheetValue(row, 'sample')
  def r1     = samplesheetValue(row, 'r1')
  def r2     = samplesheetValue(row, 'r2')
  def cram   = samplesheetValue(row, 'cram')
  def crai   = samplesheetValue(row, 'crai')
  def bam    = samplesheetValue(row, 'bam')
  def bai    = samplesheetValue(row, 'bai')
  def m2_vcf = samplesheetValue(row, 'm2_vcf')
  def ds_vcf = samplesheetValue(row, 'ds_vcf')
  def ct_vcf = samplesheetValue(row, 'ct_vcf')

  if (!sample) {
    failSamplesheet(rowNumber, "missing required 'sample' value")
  }
  if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
    failSamplesheet(rowNumber, "sample '${sample}' contains unsupported characters; use letters, numbers, '.', '_' or '-'")
  }
  if (!seenSamples.add(sample)) {
    failSamplesheet(rowNumber, "duplicate sample '${sample}'")
  }

  def modes = []
  if (r1 || r2) modes << 'fastq'
  if (cram || crai) modes << 'cram'
  if (bam || bai) modes << 'bam'
  if (m2_vcf || ds_vcf || ct_vcf) modes << 'vcf'
  if (modes.size() != 1) {
    failSamplesheet(rowNumber, "sample '${sample}' must define exactly one input mode; found ${modes ?: 'none'}")
  }

  switch (modes[0]) {
    case 'fastq':
      requireSamplesheetPath(rowNumber, sample, 'r1', r1)
      requireSamplesheetPath(rowNumber, sample, 'r2', r2)
      return ['fastq', sample, file(r1), file(r2)]
    case 'cram':
      requireSamplesheetPath(rowNumber, sample, 'cram', cram)
      requireSamplesheetPath(rowNumber, sample, 'crai', crai)
      return ['cram', sample, file(cram), file(crai)]
    case 'bam':
      requireSamplesheetPath(rowNumber, sample, 'bam', bam)
      requireSamplesheetPath(rowNumber, sample, 'bai', bai)
      return ['cram', sample, file(bam), file(bai)]
    case 'vcf':
      requireVcfIndexIfCompressed(rowNumber, sample, 'm2_vcf', m2_vcf)
      requireVcfIndexIfCompressed(rowNumber, sample, 'ds_vcf', ds_vcf)
      requireVcfIndexIfCompressed(rowNumber, sample, 'ct_vcf', ct_vcf)
      return ['vcf', sample, file(m2_vcf), file(ds_vcf), file(ct_vcf)]
    default:
      failSamplesheet(rowNumber, "unsupported input mode '${modes[0]}' for sample '${sample}'")
  }
}

// QC PROCESSES
process RAW_FASTQC {
  label 'cpu_low'
  publishDir "${params.outdir}/qc/01_raw_qc", mode: 'copy'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip

  script:
  """
  module load fastqc
  fastqc -t ${task.cpus} -q ${r1} ${r2}
  """
}

process FASTP {
  label 'cpu_med'
  publishDir "${params.outdir}/qc/02_trimmed_qc", mode: 'copy', pattern: '*.fastp.*'
  publishDir "${params.outdir}/fastq_trim", mode: 'copy', pattern: '*.trim.fq.gz'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    tuple val(sample), path("${sample}_R1.trim.fq.gz"), path("${sample}_R2.trim.fq.gz"), emit: trimmed
    path "${sample}.fastp.html", emit: html
    path "${sample}.fastp.json", emit: json

  script:
  """
  module load fastp
  fastp -w ${task.cpus} \
    -i ${r1} -I ${r2} \
    -o ${sample}_R1.trim.fq.gz -O ${sample}_R2.trim.fq.gz \
    --trim_poly_g --low_complexity_filter --overrepresentation_analysis \
    -h ${sample}.fastp.html \
    -j ${sample}.fastp.json
  """
}

process TRIMMED_FASTQC {
  label 'cpu_low'
  publishDir "${params.outdir}/qc/02_trimmed_qc", mode: 'copy'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip

  script:
  """
  module load fastqc
  fastqc -t ${task.cpus} -q ${r1} ${r2}
  """
}

process KRAKEN2 {
  label 'cpu_med'
  publishDir "${params.outdir}/qc/03_kraken_reports", mode: 'copy'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    path "${sample}.kreport", emit: kreport
    path "${sample}.out", emit: out

  script:
  """
  # Launching from an existing conda/micromamba env can trigger deactivate hooks
  # that reference unset vars. Keep nounset off through activation.
  set +u
  source ~/.bashrc
  micromamba activate ${params.kraken_env}
  set -u
  kraken2 --db ${params.kraken_db} --threads ${task.cpus} --paired --gzip-compressed \
    --output ${sample}.out \
    --report ${sample}.kreport \
    ${r1} ${r2}
  """
}

process MULTIQC {
  label 'cpu_low'
  publishDir "${params.outdir}/qc/04_final_report", mode: 'copy'

  input:
    path '*' // Collect all possible inputs dynamically into current working dir

  output:
    path "Exome_Master_QC_Report.html"
    path "Exome_Master_QC_Report_data"

  script:
  """
  module load multiqc
  multiqc . -o . \
    --filename "Exome_Master_QC_Report.html" \
    --title "Exome Project: Quality & Contamination Summary"
  """
}

// PIPELINE PROCESSES
process FQ2BAM {
  label 'gpu'
  container params.pb_container
  publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.markdup.cram*'

  input:
    tuple val(sample), path(r1), path(r2)

  output:
    tuple val(sample), path("${sample}.markdup.cram"), path("${sample}.markdup.cram.crai"), path("${sample}.recal.txt")

  script:
  def bwaOptions = ((params.bwa_opts ?: '-Y') + ' -K 10000000').trim()
  """
  pbrun fq2bam \
    --num-gpus ${params.num_gpus} \
    --ref ${params.ref} \
    --in-fq ${r1} ${r2} \
    --read-group-id ${sample}.rg1 \
    --read-group-sm ${sample} \
    --read-group-lb ${params.rg_lb} \
    --read-group-pl ${params.rg_pl} \
    --bwa-options "${bwaOptions}" \
    --out-bam ${sample}.markdup.cram \
    --out-recal-file ${sample}.recal.txt \
    --knownSites ${params.mills} \
    --knownSites ${params.dbsnp} \
    --gpuwrite \
    --gpusort
  """
}

process APPLY_BQSR_GPU {
  label 'gpu'
  container params.pb_container
  publishDir "${params.outdir}/cram_bqsr", mode: 'copy', pattern: '*.bqsr.cram*'

  input:
    tuple val(sample), path(cram), path(crai), path(recal)

  output:
    tuple val(sample), path("${sample}.bqsr.cram"), path("${sample}.bqsr.cram.crai")

  script:
  """
  pbrun applybqsr \
    --num-gpus ${params.num_gpus} \
    --ref ${params.ref} \
    --in-bam ${cram} \
    --in-recal-file ${recal} \
    --out-bam ${sample}.bqsr.cram
  """
}

process ALIGN_STATS {
  label 'cpu_med'
  publishDir "${params.outdir}/qc_stats", mode: 'copy'

  input:
    tuple val(sample), path(cram), path(crai)

  output:
    path "${sample}.flagstat.txt", emit: flagstat
    path "${sample}.idxstats.txt", emit: idxstats
    path "${sample}.samtools_depth.txt", emit: depth
    path "${sample}.mosdepth*", emit: mosdepth
    path "${sample}.bedtools_coverage.txt", emit: bedtools

  script:
  def intArg = params.exome_bed ? "-b ${params.exome_bed}" : ''
  """
  set -euo pipefail
  module load samtools bedtools
  
  samtools flagstat -@ ${task.cpus} ${cram} > ${sample}.flagstat.txt
  samtools idxstats ${cram} > ${sample}.idxstats.txt
  samtools depth ${intArg} ${cram} > ${sample}.samtools_depth.txt

  # mosdepth is provided via micromamba on this cluster.
  set +u
  source ~/.bashrc
  micromamba activate mosdepth
  set -u
  module load samtools bedtools
  export REF_PATH=${params.ref}
  # MultiQC uses the summary/dist outputs; per-base output is expensive and not
  # needed for routine QC here.
  mosdepth --no-per-base -t ${task.cpus} -f ${params.ref} ${intArg} ${sample}.mosdepth ${cram}

  if [[ -n "${params.exome_bed}" ]]; then
    ref_fai="${params.ref}.fai"
    if [[ ! -f "\$ref_fai" ]]; then
      echo "Reference FASTA index not found for bedtools sorted mode: \$ref_fai" >&2
      exit 1
    fi

    export CRAM_REFERENCE=${params.ref}
    bedtools sort -faidx "\$ref_fai" -i ${params.exome_bed} > ${sample}.targets.sorted.bed

    # bedtools coverage has no thread flag. We prefilter/decode the CRAM with
    # threaded samtools view, then feed sorted alignments to the memory-efficient
    # bedtools sweep algorithm.
    samtools view -u -@ ${task.cpus} -T ${params.ref} -M -L ${sample}.targets.sorted.bed ${cram} \
    | bedtools coverage -sorted -g "\$ref_fai" -a ${sample}.targets.sorted.bed -b stdin \
      > ${sample}.bedtools_coverage.txt
  else
    echo 'No external bed provided for bedtools' > ${sample}.bedtools_coverage.txt
  fi
  """
}

process TARGET_COVERAGE_METRICS {
  label 'cpu_low'
  publishDir "${params.outdir}/qc_target_coverage", mode: 'copy'

  input:
    tuple val(sample), path(cram), path(crai)

  output:
    path "${sample}.samtools_target_coverage.txt", emit: samtools_cov, optional: true
    path "${sample}.samtools_bedcov.raw.txt", emit: bedcov_raw, optional: true
    path "${sample}.samtools_bedcov.mean_depth.tsv", emit: bedcov_mean, optional: true
    path "${sample}.bait.interval_list"
    path "${sample}.target.interval_list"
    path "${sample}.hs_metrics.txt", emit: hs_metrics
    path "${sample}.hs_per_target_coverage.txt", emit: hs_per_target, optional: true

  script:
  def emitAuxMetrics = params.target_cov_emit_aux ? 'true' : 'false'
  def perTargetArg = params.target_cov_per_target ? "--PER_TARGET_COVERAGE ${sample}.hs_per_target_coverage.txt" : ''
  """
  set -euo pipefail
  set +u
  source ~/.bashrc
  set -u
  module load samtools apptainer

  ref="${params.ref}"
  case "\$ref" in
    *.fasta) ref_dict="\${ref%.fasta}.dict" ;;
    *.fa)    ref_dict="\${ref%.fa}.dict" ;;
    *)       ref_dict="\${ref}.dict" ;;
  esac

  if [[ ! -f "\$ref_dict" ]]; then
    echo "Reference dictionary not found: \$ref_dict" >&2
    exit 1
  fi

  if [[ "${emitAuxMetrics}" == "true" ]]; then
    samtools depth -a -b ${params.exome_bed} ${cram} \
    | awk 'BEGIN{OFS="\\t"; total=0; sum=0; ge1=0; ge10=0; ge20=0; ge30=0; ge50=0} \
         {depth=\$3; total++; sum+=depth; if(depth>=1) ge1++; if(depth>=10) ge10++; if(depth>=20) ge20++; if(depth>=30) ge30++; if(depth>=50) ge50++} \
         END{ \
           print "metric","value"; \
           print "bases_in_targets",total; \
           print "mean_depth",(total>0 ? sum/total : 0); \
           print "pct_bases_ge_1x",(total>0 ? 100*ge1/total : 0); \
           print "pct_bases_ge_10x",(total>0 ? 100*ge10/total : 0); \
           print "pct_bases_ge_20x",(total>0 ? 100*ge20/total : 0); \
           print "pct_bases_ge_30x",(total>0 ? 100*ge30/total : 0); \
           print "pct_bases_ge_50x",(total>0 ? 100*ge50/total : 0); \
         }' > ${sample}.samtools_target_coverage.txt
    samtools bedcov ${params.exome_bed} ${cram} > ${sample}.samtools_bedcov.raw.txt

    awk 'BEGIN{OFS="\\t"; print "chrom","start","end","name","bases","total_depth","mean_depth"} \
         {name=(NF>=4 ? \$4 : "."); bases=\$3-\$2; total=\$NF; mean=(bases>0 ? total/bases : 0); print \$1,\$2,\$3,name,bases,total,mean}' \
        ${sample}.samtools_bedcov.raw.txt > ${sample}.samtools_bedcov.mean_depth.tsv
  fi

  apptainer exec \
    --cleanenv \
    ${params.apptainer_bind_args} \
    ${params.gatk_container} \
    gatk BedToIntervalList \
    -I ${params.exome_bait_bed} \
    -O ${sample}.bait.interval_list \
    -SD "\$ref_dict"

  apptainer exec \
    --cleanenv \
    ${params.apptainer_bind_args} \
    ${params.gatk_container} \
    gatk BedToIntervalList \
    -I ${params.exome_bed} \
    -O ${sample}.target.interval_list \
    -SD "\$ref_dict"

  apptainer exec \
    --cleanenv \
    ${params.apptainer_bind_args} \
    ${params.gatk_container} \
    gatk CollectHsMetrics \
    -I ${cram} \
    -O ${sample}.hs_metrics.txt \
    -R ${params.ref} \
    --BAIT_INTERVALS ${sample}.bait.interval_list \
    --TARGET_INTERVALS ${sample}.target.interval_list ${perTargetArg}
  """
}

process MUTECTCALLER {
  label 'gpu'
  container params.pb_container
  publishDir "${params.outdir}/vcf_raw", mode: 'copy', pattern: '*.unfiltered.vcf.gz*'

  input:
    val pon_vcf
    tuple val(sample), path(bqsr_cram), path(bqsr_crai)

  output:
    tuple val(sample), path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), path("${sample}.f1r2.tar.gz"), path("${sample}.unfiltered.vcf.gz.stats")

  script:
  def intArg = params.exome_bed ? "--interval-file ${params.exome_bed}" : ''
  def ponArg = (pon_vcf && pon_vcf != '') ? "--pon ${pon_vcf}" : ''
  def extraArgs = params.mutect_extra_args ?: ''
  def sensitivityArgs = resolvedMutectLowLod ? '--initial-tumor-lod 1.0 --tumor-lod-to-emit 1.0' : ''
  """
  pbrun mutectcaller \
    --num-gpus ${params.num_gpus} \
    --ref ${params.ref} \
    --tumor-name ${sample} \
    --in-tumor-bam ${bqsr_cram} \
    --mutect-germline-resource ${params.mutect_germline_resource} \
    --out-vcf ${sample}.unfiltered.vcf.gz \
    --mutect-f1r2-tar-gz ${sample}.f1r2.tar.gz \
    ${intArg} \
    ${ponArg} \
    ${sensitivityArgs} \
    ${extraArgs}
  """
}

process GET_PILEUPS {
  label 'cpu_low'
  container params.gatk_container
  publishDir "${params.outdir}/vcf_filtered", mode: 'copy', pattern: '*.pileups.table*'
  input:
    tuple val(sample), path(cram), path(crai)

  output:
    tuple val(sample), path("${sample}.pileups.table")

  script:
  def intArg = params.exome_bed ? "-L ${params.exome_bed}" : ''
  """
  gatk GetPileupSummaries \
    -R ${params.ref} \
    -I ${cram} \
    -V ${params.contamination_sites_vcf} \
    ${intArg} \
    -O ${sample}.pileups.table
  """
}

process CALC_CONTAM {
  label 'cpu_low'
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
  label 'cpu_low'
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
  label 'cpu_low'
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
    --stats ${stats} \
    --contamination-table ${contam} \
    --ob-priors ${artifacts} \
    -O ${sample}.filtered.vcf.gz
  """
}

process STANDARDIZE_CALLER_VCF {
  label 'cpu_low'
  publishDir "${params.outdir}/vcf_standardized", mode: 'copy', pattern: '*.std.vcf.gz*'

  input:
    tuple val(sample), val(caller_id), path(vcf)

  output:
    tuple val(sample), val(caller_id), path("${sample}.${caller_id}.std.vcf.gz"), path("${sample}.${caller_id}.std.vcf.gz.tbi")

  script:
  """
  set -euo pipefail
  module load bcftools

  sample_count=\$(bcftools query -l ${vcf} | wc -l | awk '{print \$1}')
  if [[ "\$sample_count" -gt 1 ]]; then
    echo "Expected single-sample tumor-only VCF, found \$sample_count samples in ${vcf}" >&2
    exit 1
  fi

  input_vcf="${vcf}"
  if [[ "\$sample_count" -eq 1 ]]; then
    printf '%s\n' "${sample}" > sample_name.txt
    bcftools reheader -s sample_name.txt -o ${sample}.${caller_id}.reheader.vcf.gz ${vcf}
    bcftools index --threads ${task.cpus} -t ${sample}.${caller_id}.reheader.vcf.gz
    input_vcf="${sample}.${caller_id}.reheader.vcf.gz"
  fi

  bcftools norm \
    --threads ${task.cpus} \
    -f ${params.ref} \
    -m -both \
    -Oz -o ${sample}.${caller_id}.std.vcf.gz \
    "\$input_vcf"

  bcftools index --threads ${task.cpus} -t ${sample}.${caller_id}.std.vcf.gz
  """
}

process HARMONIZE_CALLER_VCF {
  label 'cpu_low'
  publishDir "${params.outdir}/vcf_harmonized", mode: 'copy', pattern: '*.harmonized.vcf.gz*'
  publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.harmonized_filter.tsv'

  input:
    tuple val(sample), val(caller_id), path(vcf), path(vcf_tbi)

  output:
    tuple val(sample), val(caller_id), path("${sample}.${caller_id}.harmonized.vcf.gz"), path("${sample}.${caller_id}.harmonized.vcf.gz.tbi"), emit: vcf
    path "${sample}.${caller_id}.harmonized_filter.tsv", emit: summary

  script:
  """
  set -euo pipefail
  module load bcftools

  ${params.python_bin} ${projectDir}/scripts/harmonize_filter_vcf.py \
    --sample ${sample} \
    --caller ${caller_id} \
    --vcf ${vcf} \
    --out ${sample}.${caller_id}.harmonized.vcf \
    --summary ${sample}.${caller_id}.harmonized_filter.tsv \
    --snv-min-af ${resolvedPostfilterSnvMinAf} \
    --indel-min-af ${resolvedPostfilterIndelMinAf} \
    --snv-min-alt-reads ${resolvedPostfilterSnvMinAltReads} \
    --indel-min-alt-reads ${resolvedPostfilterIndelMinAltReads} \
    --snv-min-dp ${resolvedPostfilterSnvMinDp} \
    --indel-min-dp ${resolvedPostfilterIndelMinDp}

  bcftools view -Oz -o ${sample}.${caller_id}.harmonized.vcf.gz ${sample}.${caller_id}.harmonized.vcf
  bcftools index --threads ${task.cpus} -t ${sample}.${caller_id}.harmonized.vcf.gz
  """
}

process VCF_STATS {
  label 'cpu_low'
  publishDir "${params.outdir}/vcf_stats", mode: 'copy'
  
  input:
    tuple val(sample), val(caller_id), path(vcf), path(tbi)
  
  output:
    path "${vcf.baseName}.stats"

  script:
  """
  module load bcftools
  bcftools stats ${vcf} > ${vcf.baseName}.stats
  """
}

process DEEPSOMATIC {
  label 'cpu_max'
  container params.deepsomatic_container
  publishDir "${params.outdir}/vcf_raw", mode: 'copy', pattern: '*.deepsomatic.vcf.gz*'
  publishDir "${params.outdir}/logs/deepsomatic", mode: 'copy', pattern: '*.deepsomatic.log'

  input:
    tuple val(sample), path(bqsr_cram), path(bqsr_crai)

  output:
    tuple val(sample), path("${sample}.deepsomatic.vcf.gz"), path("${sample}.deepsomatic.vcf.gz.tbi")

  script:
  def bedArg = params.exome_bed ? "--regions=${params.exome_bed}" : ''
  def extraArgs = params.deepsomatic_extra_args ?: ''
  def desiredModel = params.deepsomatic_model_type ?: 'WES_TUMOR_ONLY'
  def fallbackModel = params.deepsomatic_fallback_model_type ?: 'WGS_TUMOR_ONLY'
  def customPonArg = !isBlank(params.pon_vcf) ? "--pon_filtering=${params.pon_vcf}" : ''
  def defaultPonArg = isBlank(params.pon_vcf) ? '--use_default_pon_filtering=true' : ''
  """
  set -euo pipefail
  help_text="\$(run_deepsomatic --help 2>&1 || true)"
  selected_model="${desiredModel}"

  if ! grep -q "${desiredModel}" <<< "\$help_text"; then
    if [[ "${desiredModel}" == "WES_TUMOR_ONLY" ]] && grep -q "${fallbackModel}" <<< "\$help_text"; then
      selected_model="${fallbackModel}"
      echo "WES_TUMOR_ONLY is not supported by this DeepSomatic container. Falling back to ${fallbackModel} with regions targeting." > ${sample}.deepsomatic.log
    else
      echo "DeepSomatic model type ${desiredModel} is not supported by this container." > ${sample}.deepsomatic.log
      exit 1
    fi
  else
    echo "DeepSomatic model type ${desiredModel} is supported by this container." > ${sample}.deepsomatic.log
  fi

  run_deepsomatic \
    --model_type=\$selected_model \
    --ref=${params.ref} \
    --reads_tumor=${bqsr_cram} \
    --output_vcf=${sample}.deepsomatic.vcf.gz \
    --num_shards=${task.cpus} \
    --sample_name_tumor="${sample}" \
    ${customPonArg} \
    ${defaultPonArg} \
    ${bedArg} \
    ${extraArgs} >> ${sample}.deepsomatic.log 2>&1
  """
}

process CLAIRS_TO {
  label 'cpu_max'
  publishDir "${params.outdir}/vcf_raw", mode: 'copy', pattern: '*.clairs.vcf.gz*'
  publishDir "${params.outdir}/logs/clairs_to", mode: 'copy', pattern: '*.clairs_to.log'

  input:
    tuple val(sample), path(bqsr_cram), path(bqsr_crai)

  output:
    tuple val(sample), path("${sample}.clairs.vcf.gz"), path("${sample}.clairs.vcf.gz.tbi")

  script:
  def bedArg = params.exome_bed ? "--bed_fn ${params.exome_bed}" : ""
  def clairsBind = params.apptainer_bind_args ?: ''
  def extraArgs = params.clairs_to_extra_args ?: ''
  def thresholdArgs = "--snv_min_af ${resolvedClairsToSnvMinAf} --indel_min_af ${resolvedClairsToIndelMinAf} --min_coverage ${resolvedClairsToMinCoverage}"
  def useCustomPon = [params.clairs_to_gnomad, params.clairs_to_dbsnp, params.clairs_to_pon, params.clairs_to_colorsdb].every { pathExists(it) }
  def ponArg = useCustomPon ? "--panel_of_normals \"${params.clairs_to_gnomad},${params.clairs_to_dbsnp},${params.clairs_to_pon},${params.clairs_to_colorsdb}\" --panel_of_normals_require_allele_matching \"True,True,False,False\"" : ''
  """
  set -euo pipefail
  set +u
  source ~/.bashrc
  set -u
  module load apptainer

  mkdir clairs_to_out
  apptainer exec \
    --cleanenv \
    ${clairsBind} \
    ${params.clairs_to_container} \
    run_clairs_to \
    --tumor_bam_fn ${bqsr_cram} \
    --ref_fn ${params.ref} \
    --sample_name ${sample} \
    --threads ${task.cpus} \
    --platform ${params.clairs_to_platform} \
    ${thresholdArgs} \
    ${ponArg} \
    --output_dir ./clairs_to_out \
    --conda_prefix /opt/micromamba/envs/clairs-to \
    ${bedArg} \
    ${extraArgs} \
    > ${sample}.clairs_to.log 2>&1
  
  # ClairS-TO outputs snv.vcf.gz and indel.vcf.gz separately.
  # We concatenate and sort them into a unified VCF file to match other tools.
  snv_records=0
  indel_records=0
  if [[ -f clairs_to_out/snv.vcf.gz ]]; then
    snv_records=\$(apptainer exec --cleanenv ${clairsBind} ${params.clairs_to_container} bcftools view -H clairs_to_out/snv.vcf.gz | wc -l | awk '{print \$1}')
  fi
  if [[ -f clairs_to_out/indel.vcf.gz ]]; then
    indel_records=\$(apptainer exec --cleanenv ${clairsBind} ${params.clairs_to_container} bcftools view -H clairs_to_out/indel.vcf.gz | wc -l | awk '{print \$1}')
  fi
  if [[ "\$snv_records" -eq 0 && "\$indel_records" -eq 0 ]]; then
    echo "WARNING: ClairS-TO produced no SNV/indel records for ${sample}" >&2
  fi

  if [[ -f clairs_to_out/snv.vcf.gz && -f clairs_to_out/indel.vcf.gz ]]; then
    apptainer exec \
      --cleanenv \
      ${clairsBind} \
      ${params.clairs_to_container} \
      bcftools concat -a -Ou clairs_to_out/snv.vcf.gz clairs_to_out/indel.vcf.gz \
    | apptainer exec \
        --cleanenv \
        ${clairsBind} \
        ${params.clairs_to_container} \
        bcftools sort -Oz -o ${sample}.clairs.vcf.gz
  elif [[ -f clairs_to_out/snv.vcf.gz ]]; then
    cp -f clairs_to_out/snv.vcf.gz ${sample}.clairs.vcf.gz
  elif [[ -f clairs_to_out/indel.vcf.gz ]]; then
    cp -f clairs_to_out/indel.vcf.gz ${sample}.clairs.vcf.gz
  else
    echo "ClairS-TO finished without producing snv.vcf.gz or indel.vcf.gz" >&2
    exit 1
  fi

  apptainer exec \
    --cleanenv \
    ${clairsBind} \
    ${params.clairs_to_container} \
    bcftools index --threads ${task.cpus} -t ${sample}.clairs.vcf.gz
  """
}

process CONCORDANCE_MERGE {
  label 'cpu_low'
  publishDir "${params.outdir}/vcf_concordant", mode: 'copy'
  publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.tsv'

  input:
    tuple val(sample), path(m2_vcf), path(m2_tbi), path(ds_vcf), path(ds_tbi), path(ct_vcf), path(ct_tbi)

  output:
    tuple val(sample), path("${sample}.concordant.vcf.gz"), path("${sample}.concordant.vcf.gz.tbi"), emit: concordant
    tuple val(sample), path("${sample}.all_callers.vcf.gz"), path("${sample}.all_callers.vcf.gz.tbi"), emit: all_callers
    path "${sample}.concordant.vcf.stats", emit: stats
    path "${sample}.caller_support.tsv", emit: support_tsv
    path "${sample}.all_callers_support.tsv", emit: all_support_tsv
    path "${sample}.caller_characteristics.tsv", emit: caller_characteristics
    path "${sample}.concordance_summary.tsv", emit: concordance_summary
    path "${sample}.pairwise_concordance.tsv", emit: pairwise_concordance

  script:
  """
  set -euo pipefail
  module load bcftools

  mkdir -p isec_out
  bcftools isec \
    -c none \
    -n+${resolvedMergeMinCallers} \
    -p isec_out \
    ${m2_vcf} \
    ${ds_vcf} \
    ${ct_vcf}

  ${params.python_bin} ${projectDir}/scripts/merge_caller_support.py \
    --sample ${sample} \
    --min-callers ${resolvedMergeMinCallers} \
    --candidate-dir isec_out \
    --output-vcf ${sample}.concordant.vcf \
    --support-tsv ${sample}.caller_support.tsv \
    --all-output-vcf ${sample}.all_callers.vcf \
    --all-support-tsv ${sample}.all_callers_support.tsv \
    --caller-characteristics-tsv ${sample}.caller_characteristics.tsv \
    --concordance-summary-tsv ${sample}.concordance_summary.tsv \
    --pairwise-concordance-tsv ${sample}.pairwise_concordance.tsv \
    --caller mutect2:${m2_vcf} \
    --caller deepsomatic:${ds_vcf} \
    --caller clairsto:${ct_vcf}

  bcftools sort -Oz -o ${sample}.concordant.vcf.gz ${sample}.concordant.vcf
  bcftools index --threads ${task.cpus} -t ${sample}.concordant.vcf.gz
  bcftools sort -Oz -o ${sample}.all_callers.vcf.gz ${sample}.all_callers.vcf
  bcftools index --threads ${task.cpus} -t ${sample}.all_callers.vcf.gz
  bcftools stats ${sample}.concordant.vcf.gz > ${sample}.concordant.vcf.stats
  """
}

process ANNOTATE_SNPEFF {
  label 'cpu_med'
  publishDir "${params.outdir}/vcf_annotated", mode: 'copy'

  input:
    tuple val(sample), path(vcf), path(tbi)

  output:
    tuple val(sample), path("${sample}.snpeff.vcf.gz"), path("${sample}.snpeff.vcf.gz.tbi"), emit: vcf
    path "${sample}.snpEff_summary.html", emit: summary
    path "${sample}.snpEff_genes.txt", emit: genes

  script:
  def snpeffJar = params.snpeff_jar ?: "/apps/snpEff/${params.snpeff_version}/bin/snpEff.jar"
  def extraArgs = params.snpeff_extra_args ?: ''
  def annotateArgs = params.annotation_bcftools_args ?: ''
  """
  set -euo pipefail
  set +u
  source ~/.bashrc
  module load snpEff/${params.snpeff_version} bcftools || true
  set -u
  
  java -Xmx${params.snpeff_java_mem} -jar ${snpeffJar} \
    ${params.snpeff_db} \
    ${extraArgs} \
    -stats ${sample}.snpEff_summary.html \
    ${vcf} \
  | bcftools sort -Oz -o ${sample}.snpeff.tmp.vcf.gz

  bcftools index --threads ${task.cpus} -t ${sample}.snpeff.tmp.vcf.gz

  if [[ -n "${annotateArgs}" ]]; then
    bcftools annotate ${annotateArgs} \
      -Oz -o ${sample}.snpeff.vcf.gz \
      ${sample}.snpeff.tmp.vcf.gz
  else
    mv ${sample}.snpeff.tmp.vcf.gz ${sample}.snpeff.vcf.gz
    mv ${sample}.snpeff.tmp.vcf.gz.tbi ${sample}.snpeff.vcf.gz.tbi
  fi

  if [[ ! -f ${sample}.snpeff.vcf.gz.tbi ]]; then
    bcftools index --threads ${task.cpus} -t ${sample}.snpeff.vcf.gz
  fi

  if [[ -f snpEff_genes.txt ]]; then
    mv snpEff_genes.txt ${sample}.snpEff_genes.txt
  else
    touch ${sample}.snpEff_genes.txt
  fi
  """
}

process VARIANT_REVIEW_TABLE {
  label 'cpu_low'
  publishDir "${params.outdir}/tables", mode: 'copy', pattern: '*.variant_review.tsv'

  input:
    tuple val(sample), path(vcf), path(tbi)

  output:
    path "${sample}.variant_review.tsv"

  script:
  """
  set -euo pipefail
  ${params.python_bin} ${projectDir}/scripts/variants_to_tsv.py \
    --sample ${sample} \
    --vcf ${vcf} \
    --out ${sample}.variant_review.tsv
  """
}



process COPY_REF_TO_CRAM_DIR {
  label 'cpu_low'
  container params.gatk_container
  publishDir "${params.outdir}/cram", mode: 'copy', overwrite: true

  input:
    path ref_fa
    val  _

  output:
    path "${ref_fa.getFileName()}", emit: fasta
    path "${ref_fa.getFileName()}.fai", emit: fai, optional: true
    path "*.dict", emit: dict, optional: true

  script:
  """
  set -euo pipefail
  ref_src="${params.ref}"
  ref_base=\$(basename "\$ref_src")
  tmp=".nf_tmp_\$ref_base"

  # Make a new file first, then replace the staged input path atomically.
  rm -f "\$tmp"
  ln -f "\$ref_src" "\$tmp" 2>/dev/null || cp -a --reflink=auto "\$ref_src" "\$tmp"
  mv -f "\$tmp" "\$ref_base"

  if [[ -f "\${ref_src}.fai" ]]; then
    cp -af "\${ref_src}.fai" "\${ref_base}.fai"
  fi

  case "\$ref_src" in
    *.fasta) ref_dict="\${ref_src%.fasta}.dict" ;;
    *.fa)    ref_dict="\${ref_src%.fa}.dict" ;;
    *)       ref_dict="\${ref_src}.dict" ;;
  esac
  if [[ -f "\$ref_dict" ]]; then
    cp -af "\$ref_dict" .
  fi
  """
}

// ---------------------- TOP-LEVEL ----------------------
workflow {
  if( !params.samplesheet )
    exit 1, "Missing required --samplesheet or SAMPLESHEET environment variable"

  if( !file(params.samplesheet).exists() )
    exit 1, "Samplesheet not found: ${params.samplesheet}"

  input_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true, sep:',')
    .collect()
    .flatMap { rows ->
      if (!rows) {
        exit 1, "Samplesheet has no data rows: ${params.samplesheet}"
      }
      def seenSamples = [] as Set
      def parsedRows = []
      rows.eachWithIndex { row, idx ->
        parsedRows << parseSamplesheetRow(row, idx + 2, seenSamples)
      }
      parsedRows
    }

  branched_input = input_ch.branch {
    fastq: it[0] == 'fastq'
    cram:  it[0] == 'cram'
    vcf:   it[0] == 'vcf'
  }

  reads_ch = branched_input.fastq.map { tuple(it[1], it[2], it[3]) }
  ready_cram_ch = branched_input.cram.map { tuple(it[1], it[2], it[3]) }
  ready_vcf_ch = branched_input.vcf.map { tuple(it[1], it[2], it[3], it[4]) }
  
  // ---------------- FASTQ PREPROCESSING ----------------
  // 1. Raw FastQC
  raw_qc = RAW_FASTQC(reads_ch)

  // 2. Fastp trimming conditionally
  if (params.skip_trimming) {
    trimmed_reads_ch = reads_ch
    fastp_json = Channel.empty()
    fastp_html = Channel.empty()
    trimmed_qc_html = Channel.empty()
    trimmed_qc_zip = Channel.empty()
  } else {
    fastp_out = FASTP(reads_ch)
    trimmed_reads_ch = fastp_out.trimmed
    fastp_json = fastp_out.json
    fastp_html = fastp_out.html

    trimmed_qc = TRIMMED_FASTQC(trimmed_reads_ch)
    trimmed_qc_html = trimmed_qc.html
    trimmed_qc_zip = trimmed_qc.zip
  }
  
  // 4. Kraken2
  k2_out = KRAKEN2(trimmed_reads_ch)

  mutect_pon_vcf_ch = Channel.value( params.mutect_pon_vcf ?: '' )
  ref_ch = Channel.of(params.ref)
  
  // Align trimmed reads
  fq2   = FQ2BAM(trimmed_reads_ch)
  fq2_done_ch = fq2.collect().map { true }
  COPY_REF_TO_CRAM_DIR(ref_ch, fq2_done_ch)
  
  // Apply BQSR
  bqsr_cram_from_fq = APPLY_BQSR_GPU(fq2)
  
  // Combine newly aligned bqsr_crams with direct BAM/CRAM inputs
  final_cram_ch = bqsr_cram_from_fq.mix(ready_cram_ch)
  
  // ---------------- METRICS & CALLING ----------------
  // Map stats
  align_stats_out = ALIGN_STATS(final_cram_ch)
  target_cov_out = TARGET_COVERAGE_METRICS(final_cram_ch)
  
  // Callers
  m2  = MUTECTCALLER(mutect_pon_vcf_ch, final_cram_ch)
  ds  = DEEPSOMATIC(final_cram_ch)
  ct  = CLAIRS_TO(final_cram_ch)

  // Mutect2 post-processing & filtering
  m2_vcf4f = m2
  pile     = GET_PILEUPS(final_cram_ch)
  contam   = CALC_CONTAM(pile)
  orient   = LEARN_ORIENT(m2_vcf4f)

  joined_m2 = orient
    .join(contam, by: 0)
    .map { row ->
      def (s, artifacts, vcf_in, vcf_tbi, stats, cont) = row
      tuple(s, artifacts, vcf_in, vcf_tbi, stats, cont)
    }

  filt_m2  = FILTER_M2(joined_m2)
  
  // Standardize and harmonize all caller VCFs after caller-native generation.
  // Mutect2 still runs FilterMutectCalls, but final comparison is driven by the
  // shared AF/ALT-read/DP filter rather than caller-native PASS alone.
  m2_for_std = filt_m2.map { s, vcf, tbi -> tuple(s, 'mutect2', vcf) }
  ds_for_std = ds.map { s, vcf, tbi -> tuple(s, 'deepsomatic', vcf) }
  ct_for_std = ct.map { s, vcf, tbi -> tuple(s, 'clairsto', vcf) }

  user_m2_for_std = ready_vcf_ch.map { s, m2, ds, ct -> tuple(s, 'mutect2', m2) }
  user_ds_for_std = ready_vcf_ch.map { s, m2, ds, ct -> tuple(s, 'deepsomatic', ds) }
  user_ct_for_std = ready_vcf_ch.map { s, m2, ds, ct -> tuple(s, 'clairsto', ct) }

  std_inputs = m2_for_std
    .mix(ds_for_std)
    .mix(ct_for_std)
    .mix(user_m2_for_std)
    .mix(user_ds_for_std)
    .mix(user_ct_for_std)

  standardized_all = STANDARDIZE_CALLER_VCF(std_inputs)
  harmonized_out = HARMONIZE_CALLER_VCF(standardized_all)
  harmonized_stats = VCF_STATS(harmonized_out.vcf)

  m2_norm = harmonized_out.vcf
    .filter { s, caller, vcf, tbi -> caller == 'mutect2' }
    .map { s, caller, vcf, tbi -> tuple(s, vcf, tbi) }

  ds_norm = harmonized_out.vcf
    .filter { s, caller, vcf, tbi -> caller == 'deepsomatic' }
    .map { s, caller, vcf, tbi -> tuple(s, vcf, tbi) }

  ct_norm = harmonized_out.vcf
    .filter { s, caller, vcf, tbi -> caller == 'clairsto' }
    .map { s, caller, vcf, tbi -> tuple(s, vcf, tbi) }
  
  // Combine standardized VCFs by sample for consensus calling.
  called_vcf_ch = m2_norm.join(ds_norm, by: 0).join(ct_norm, by: 0)

  // ---------------- CONCORDANCE & ANNOTATION ----------------
  concord = CONCORDANCE_MERGE(called_vcf_ch)
  
  snpeff_out = ANNOTATE_SNPEFF(concord.concordant)
  review_tables = VARIANT_REVIEW_TABLE(snpeff_out.vcf)

  // ---------------- MULTIQC ----------------
  all_qc_files = raw_qc.html.mix(raw_qc.zip)
    .mix(trimmed_qc_html).mix(trimmed_qc_zip)
    .mix(fastp_json).mix(fastp_html)
    .mix(k2_out.kreport).mix(k2_out.out)
    .mix(align_stats_out.flagstat).mix(align_stats_out.idxstats)
    .mix(align_stats_out.mosdepth).mix(align_stats_out.bedtools)
    .mix(target_cov_out.hs_metrics)
    .mix(harmonized_stats).mix(concord.stats)
    .collect()
    
  MULTIQC(all_qc_files)

  emit:
    annotated_vcf = snpeff_out.vcf
    review_tsv = review_tables
    concordant_vcf = concord.concordant
}
