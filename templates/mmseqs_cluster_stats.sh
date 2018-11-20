mmseqs createtsv \
  "${seq}/db" \
  "${seq}/db" \
  "${clusters}/db" \
  "${clusters}.tsv" \
  --threads ${task.cpus}

sed -i '1i cluster\tmember' "${clusters}.tsv"


mmseqs result2repseq \
  "${seq}/db" \
  "${clusters}/db" \
  "${clusters}_rep" \
  --threads ${task.cpus}

mmseqs result2flat \
  "${seq}/db" \
  "${seq}/db" \
  "${clusters}_rep" \
  "${clusters}_rep.fasta" \
  --use-fasta-header


# Do an all vs centroid alignment for each cluster.
mmseqs align \
  "${seq}/db" \
  "${seq}/db" \
  "${clusters}/db" \
  align \
  -a \
  --threads ${task.cpus}

mmseqs convertalis \
  "${seq}/db" \
  "${seq}/db" \
  align \
  "${clusters}_stats.tsv" \
  --threads ${task.cpus} \
  --format-mode 0 \
  --format-output "query target evalue qcov tcov gapopen pident nident mismatch raw bits qstart qend tstart tend qlen tlen alnlen"

sed -i '1i query\ttarget\tevalue\tqcov\ttcov\tgapopen\tpident\tnident\tmismatch\traw\tbits\tqstart\tqend\ttstart\ttend\tqlen\ttlen\talnlen' "${clusters}_stats.tsv"
