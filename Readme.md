# Parabricks Mutect2 (Tumor‑only) — Nextflow Pipeline

This repository runs **tumor‑only somatic variant calling** using NVIDIA **Parabricks** (GPU‑accelerated) and **GATK** helper tools. It aligns FASTQs → CRAM with `pbrun fq2bam`, then calls variants with `pbrun mutectcaller`, and finally applies **GetPileupSummaries → CalculateContamination → LearnReadOrientationModel → FilterMutectCalls → SelectVariants** to produce a *PASS‑only* VCF.

---

## Pipeline at a glance

**Inputs (per sample)**

* Paired‑end FASTQs (`r1`, `r2`)
* Reference genome (FASTA + index & dict)
* Optional:

  * Exome BED intervals (`params.exome_bed`)
  * Panel of Normals VCF (`params.pon_vcf`)
  * Germline resource (gnomAD AF; `params.gnomad_af`)

**Core steps / processes**

1. **FQ2BAM** *(GPU, Parabricks)* — `pbrun fq2bam`

   * Aligns with BWA‑MEM (`-Y -K 10,000,000`), marks dups, BQSR.
   * Outputs: `{sample}.markdup.cram`, `.crai`, `.recal.txt`.
2. **MUTECTCALLER** *(GPU, Parabricks)* — `pbrun mutectcaller`

   * Uses `--mutect-germline-resource` (= gnomAD AF).
   * Optional: `--interval-file` (exome BED), `--pon` (PON VCF).
   * Outputs: `{sample}.unfiltered.vcf.gz`, `{sample}.f1r2.tar.gz`, stats.
3. **GET_PILEUPS** *(CPU, GATK)* — `gatk GetPileupSummaries`
4. **CALC_CONTAM** *(CPU, GATK)* — `gatk CalculateContamination`
5. **LEARN_ORIENT** *(CPU, GATK)* — `gatk LearnReadOrientationModel`
6. **FILTER_M2** *(CPU, GATK)* — `gatk FilterMutectCalls`
7. **SELECT_PASS** *(CPU, GATK)* — `gatk SelectVariants --exclude-filtered`
8. **COPY_REF_TO_CRAM_DIR** *(CPU)* — Copies/creates `.fai` and `.dict` into the CRAM publish directory for convenience.

**Outputs**

```
{outdir}/
├─ cram/
│   ├─ {sample}.markdup.cram
│   ├─ {sample}.markdup.cram.crai
│   ├─ {ref}.{fai,dict}         # created/copied by helper step
├─ vcf/
│   ├─ {sample}.unfiltered.vcf.gz{,.tbi,.stats}
│   ├─ {sample}.filtered.vcf.gz{,.tbi}
│   └─ {sample}.pass.vcf.gz{,.tbi}
```

---

## Input samplesheet (CSV)

The workflow **expects a CSV** (comma‑separated) with header **`sample,r1,r2`**. Example (`sample_2.csv`):

```csv
sample,r1,r2
sample_1,/path/to/fastq/sample.R1.fq.gz,/path/to/ALL_FF67/fastq/sample_R2.trim.fq.gz
```

> If a cell is empty or the header doesn’t exactly match, the pipeline will **stop** with a clear error.

Pass the file path via `--samplesheet`.

---

## Parameters (from `main.nf` / `nextflow.config`)

| Param                         | Meaning                         | Notes                                       |
| ----------------------------- | ------------------------------- | ------------------------------------------- |
| `--samplesheet`               | CSV with columns `sample,r1,r2` | required                                    |
| `--outdir`                    | Output base directory           | default `results`                           |
| `params.ref`                  | Reference FASTA                 | required; must match all known sites        |
| `params.mills`                | Mills indels VCF                | for BQSR                                    |
| `params.dbsnp`                | dbSNP VCF                       | for BQSR                                    |
| `params.gnomad_af`            | gnomAD AF VCF                   | Mutect2 germline resource                   |
| `params.exome_bed`            | BED intervals                   | optional; used in `fq2bam` and Mutect steps |
| `params.pon_vcf`              | PON VCF                         | optional; Mutect `--pon`                    |
| `params.num_gpus`             | GPUs per Parabricks step        | default `2` in `main.nf`                    |
| `params.rg_lb`/`params.rg_pl` | Read group fields               | e.g., `sureselect_v8` / `ILLUMINA`          |
| `params.pb_container`         | Parabricks container path       | set in `nextflow.config`                    |
| `params.gatk_container`       | GATK container path             | set in `nextflow.config`                    |

**GPU vs CPU**

* `FQ2BAM` and `MUTECTCALLER` are labeled `gpu`.
* The rest run on CPU (`gatk` container).

---

## How to run

**Slurm/GPU example (from `mutact.sh`, generalized):**

```bash
module load apptainer
module load nextflow
export NXF_SINGULARITY_CMD=apptainer

nextflow run . -profile slurm -resume \
  --samplesheet path/to/sample.csv \
  --outdir results/parabricks_nf_out \
  --num_gpus 1 \
  --rg_lb sureselect_v8 --rg_pl ILLUMINA
```

> In your original script the path to `main.nf` and the samplesheet are **absolute**. For a public repo, prefer **relative paths** (as above) and move cluster specifics into a **private config**.

---

## Reproducibility notes

* Pin exact **container** versions (e.g. Parabricks `4.5.1`, GATK `4.6.1.0`).
* Keep the same **reference build** across FASTA and all VCFs (hg38 in your example).
* Save the effective `nextflow.config` (`-with-trace`, `-with-report`, `-with-timeline`) per run.

---

## License

Add a license (MIT/Apache‑2.0) to clarify reuse.

