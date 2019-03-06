# Search the profiles against the profile consensus sequences.
# Uses an iterative strategy.

# Params
# INDB
# NCPUS
# OUTDB

mkdir -p tmp
mkdir search_result

mmseqs search \
  "${INDB}/db" \
  "${INDB}/db_consensus" \
  search_result/db \
  tmp \
  --threads ${NCPUS} \
  --max-seqs 300 \
  -c 0.8 \
  --cov-mode 1 \
  --start-sens 5 \
  --sens-steps 2 \
  -s 7.0 \
  --add-self-matches \
  --num-iterations 2

# Cluster the matches of profiles vs consensus sequences.
mkdir -p ${OUTDB}
mmseqs clust \
  "${INDB}/db" \
  search_result/db \
  "${OUTDB}/db" \
  --threads ${NCPUS} \
  --cluster-mode 2

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
