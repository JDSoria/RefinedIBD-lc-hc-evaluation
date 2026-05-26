#!/bin/bash
#SBATCH --job-name=IBD_HC
#SBATCH --partition=mono
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --output=logHC/IBD_%A_%a.out
#SBATCH --error=logHC/IBD_%A_%a.err

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:-${1:?Usage: bash RunRefinedHC.sh <chr>}}

BASE=${BASE:-/home/dsoria/PoblAR/LocalAncestryMuestrasPiloto/Script/Input}
VCF=${VCF:-/home/dsoria/PoblAR/RefinedIBD/HC_clean_AF01/chr${CHR}_clean.vcf.gz}
MAP=${MAP:-${BASE}/MapasGeneticos/chr${CHR}.b38.gmap.plink}
OUTDIR=${OUTDIR:-/home/dsoria/PoblAR/RefinedIBD/HC}
IMAGE=${IMAGE:-${BASE}/RefinedIBD/refined_ibd.sif}

mkdir -p "${OUTDIR}" logHC

[[ -f "${VCF}" ]] || { echo "Missing VCF: ${VCF}" >&2; exit 1; }
[[ -f "${MAP}" ]] || { echo "Missing genetic map: ${MAP}" >&2; exit 1; }
[[ -f "${IMAGE}" ]] || { echo "Missing Singularity image: ${IMAGE}" >&2; exit 1; }

echo "Running RefinedIBD HC chromosome ${CHR}"

singularity exec --bind /home/dsoria:/home/dsoria \
  "${IMAGE}" \
  java -Xmx38g -jar /opt/refined-ibd.17Jan20.102.jar \
  gt="${VCF}" \
  map="${MAP}" \
  out="${OUTDIR}/chr${CHR}_ibd" \
  window=20 \
  trim=0.3

echo "RefinedIBD chr${CHR} DONE (HC)"

IBD_GZ="${OUTDIR}/chr${CHR}_ibd.ibd.gz"
MERGED_OUT="${OUTDIR}/chr${CHR}_ibd_merged.ibd"
TMP="${OUTDIR}/chr${CHR}_ibd.tmp"

[[ -f "${IBD_GZ}" ]] || { echo "Missing RefinedIBD output: ${IBD_GZ}" >&2; exit 1; }

trap 'rm -f "${TMP}"' EXIT

echo "Merging segments chr${CHR} (HC)"

gunzip -c "${IBD_GZ}" > "${TMP}"

# merge-ibd-segments reads the RefinedIBD IBD file from stdin.
singularity exec --bind /home/dsoria:/home/dsoria \
  "${IMAGE}" \
  java -jar /opt/merge-ibd-segments.17Jan20.102.jar \
  "${VCF}" "${MAP}" 0.01 1 \
  < "${TMP}" \
  > "${MERGED_OUT}"

rm -f "${TMP}"
trap - EXIT

echo "Chromosome ${CHR} DONE (HC)"
echo "Merged output: ${MERGED_OUT}"
