# Params

# QUERY
# TARGET
# RESULTS
# OUTDB

mkdir -p "${OUTDB}"

mmseqs result2profile \
  "${QUERY}/db" \
  "${TARGET}/db" \
  "${RESULTS}/db" \
  "${OUTDB}/db" \
  --threads ${NCPUS}

#mkdir -p tmp
#mmseqs createindex "${OUTDB}/db_consensus" tmp --threads "${NCPUS}"

#rm -rf -- tmp
