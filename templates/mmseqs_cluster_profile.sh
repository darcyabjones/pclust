# Search the profiles against the profile consensus sequences.
# Uses an iterative strategy.

# Params
# INDB
# NCPUS
# OUTDB

mkdir -p tmp
mkdir search_result

# mmseqs touchdb "${INDB}/db_consensus" --threads "${NCPUS}"

mmseqs search \
  "${INDB}/db" \
  "${INDB}/db_consensus" \
  search_result/db \
  tmp \
  --threads ${NCPUS} \
  --max-seqs 300 \
  -c 0.8 \
  --cov-mode 0 \
  --start-sens 4 \
  --sens-steps 3 \
  -s 7.5 \
  -e 0.00001 \
  --e-profile 0.01 \
  --add-self-matches \
  --db-load-mode 3

# Cluster the matches of profiles vs consensus sequences.
mkdir -p ${OUTDB}
mmseqs clust \
  "${INDB}/db" \
  search_result/db \
  "${OUTDB}/db" \
  --threads ${NCPUS} \
  --cluster-mode 0

mmseqs convertalis \
  ${INDB}/db \
  ${INDB}/db_consensus \
  search_result/db \
  "${OUTDB}_matches.tsv" \
  --threads ${NCPUS} \
  --format-mode 0 \
  --format-output "query,target,evalue,qcov,tcov,gapopen,pident,nident,mismatch,raw,bits,qstart,qend,tstart,tend,qlen,tlen,alnlen"

sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "${OUTDB}_matches.tsv"

rm -rf -- tmp
rm -rf -- search_result
