mmseqs createtsv \
  "${seq}/db" \
  "${seq}/db" \
  "${clusters}/db" \
  "${clusters}.tsv"

sed -i '1i cluster\tmember' "${clusters}.tsv"


mmseqs result2repseq \
  "${seq}/db" \
  "${clusters}/db" \
  "${clusters}_rep"

mmseqs result2flat \
  "${seq}/db" \
  "${seq}/db" \
  "${clusters}_rep" \
  "${clusters}_rep.fasta" \
  --use-fasta-header


# Do an all vs centroid alignment for each cluster.
mmseqs align "${seq}/db" "${seq}/db" "${clusters}/db" align -a
mmseqs convertalis "${seq}/db" "${seq}/db" align "${clusters}_stats.tsv"
sed -i '1i query\ttarget\tident\tlength\tmismatch\tngap\tqstart\tqend\ttstart\ttend\tevalue\tbitscore' "${clusters}_stats.tsv"
