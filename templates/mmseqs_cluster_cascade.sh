# Params
# INDB
# OUTDB
# NCPUS

mkdir -p "${OUTDB}"
mkdir -p tmp

mmseqs cluster \
  "${INDB}/db" \
  "${CASCADE}/db" \
  tmp \
  --threads ${NCPUS} \
  --min-seq-id 0.0 \
  -c 0.7 \
  --cov-mode 0 \
  -s 7.5 \
  --cluster-steps 5 \
  --cluster-mode 0

rm -rf -- tmp
