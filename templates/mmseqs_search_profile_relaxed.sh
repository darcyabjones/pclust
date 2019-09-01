# Params

# PROFILE
# TARGET
# OUTDB
# NCPUS

mkdir -p tmp
mkdir -p "${OUTDB}"

mmseqs search \
  "${PROFILE}/db" \
  "${TARGET}/db" \
  "${OUTDB}/db" \
  tmp \
  -a \
  --threads ${NCPUS} \
  --max-seqs 300 \
  -e 0.00001 \
  --e-profile 0.01 \
  --start-sens 4.0 \
  --sens-steps 3 \
  -s 7.5 \
  --rescore-mode 1 \
  --realign \
  --db-load-mode 3

rm -rf -- tmp
