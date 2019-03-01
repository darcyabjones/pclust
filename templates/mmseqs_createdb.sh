# Params:
# FASTA
# OUTDB

mkdir -p "${OUTDB}"
mmseqs createdb "${FASTA}" "${OUTDB}/db" --max-seq-len 14000
