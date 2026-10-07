#!/usr/bin/env bash
# Why real reads are soft-clipped: every clip in the HG002 chr20 slice, one cause per clipped read.
#   REF       GRCh38 no-alt analysis set FASTA, bwa-mem2 indexed
#   SLICE     chr20:38.5-40.2 Mb slice of the GIAB HG002 NovaSeq PCR-free 35x BAM (bwa-mem2 2.2.1)
#   Q100_VCF  GIAB HG002 T2T-Q100 v1.1 benchmark VCF (GRCh38_HG2-T2TQ100-V1.1.vcf.gz)
set -euo pipefail
cd "$(dirname "$0")"
: "${REF:?}" "${SLICE:?}" "${Q100_VCF:?}"
python3 collect.py                                    # clips.pkl, clips_ge20.fa
bwa-mem2 mem -t 32 "$REF" clips_ge20.fa 2> realign.log | samtools view -b -o clips_ge20.bam -
python3 classify.py                                   # first pass (adapter check missed the R2 adapter)
python3 classify2.py                                  # cls2.pkl: adapter, site, chimera, foreign, polyG, badend, hiQerr
python3 foreign.py                                    # foreign_cls.pkl: what the foreign clips are
python3 foreign_sw.py                                 # foreign_sw.pkl: local alignment of the unexplained ones
python3 rest_fasta.py
bwa-mem2 mem -t 32 -k 13 -T 20 -B 2 "$REF" rest_ge20.fa 2> rest.log | samtools view -b -o rest_ge20.bam -
python3 rest_realign_count.py
python3 summary.py                                    # one cause per clipped read
python3 nonerr_mask.py                                # nonerr_mask.pkl, used by ../qmodel/e4m.py and e9m.py
