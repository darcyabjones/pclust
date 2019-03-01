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
  --threads ${NCPUS} \
  --max-seqs 200 \
  -e 0.00001 \
  --e-profile 0.01 \
  -c 0.2 \
  --start-sens 5.0 \
  --sens-steps 2 \
  -s 7.0 \
  --rescore-mode 1 \
  --num-iterations 2

rm -rf -- tmp
