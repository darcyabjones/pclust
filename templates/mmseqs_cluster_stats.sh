# Params:
# SEQS
# CLUSTERS
# NCPUS


mmseqs createtsv \
  "${SEQS}/db" \
  "${SEQS}/db" \
  "${CLUSTERS}/db" \
  "${CLUSTERS}.tsv" \
  --threads ${NCPUS}

sed -i '1i cluster\tmember' "${CLUSTERS}.tsv"


mmseqs result2repseq \
  "${SEQS}/db" \
  "${CLUSTERS}/db" \
  "${CLUSTERS}_rep" \
  --threads ${NCPUS}

mmseqs result2flat \
  "${SEQS}/db" \
  "${SEQS}/db" \
  "${CLUSTERS}_rep" \
  "${CLUSTERS}_rep.fasta" \
  --use-fasta-header


# Do an all vs centroid alignment for each cluster.
mmseqs align \
  "${SEQS}/db" \
  "${SEQS}/db" \
  "${CLUSTERS}/db" \
  align \
  -a \
  --threads ${NCPUS}

mmseqs convertalis \
  "${SEQS}/db" \
  "${SEQS}/db" \
  align \
  "${CLUSTERS}_stats.tsv" \
  --threads ${NCPUS} \
  --format-mode 0 \
  --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen"

sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "${CLUSTERS}_stats.tsv"
