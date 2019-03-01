mkdir -p "${outdb}"
mkdir -p "tmp"

mmseqs linclust \
  "${indb}/db" \
  "${outdb}/db" \
  tmp \
  --min-seq-id 0.90 \
  -c 0.8 \
  --cov-mode 0

rm -rf -- tmp
