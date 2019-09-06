#!/usr/bin/env bash

SEQS="$(cat /dev/stdin | tr '\0' '\n')"
NSEQS=$(echo "${SEQS}" | grep -c "^>" )

if [ ${NSEQS} -le 1 ]
then
  echo "${SEQS}"
elif [ ${NSEQS} -le 200 ]
then
  # Good but slow.
  echo "${SEQS}" \
  | mafft \
      --anysymbol \
      --inputorder \
      --quiet \
      --bl 30 \
      --genafpair \
      --maxiterate 1000 \
      --amino \
      --ep 0.0 \
      -
elif [ ${NSEQS} -le 10000 ]
then
  # OK but faster
  echo "${SEQS}" \
  | mafft \
      --anysymbol \
      --inputorder \
      --quiet \
      --bl 30 \
      --6merpair \
      --retree 2 \
      --maxiterate 1000 \
      --amino \
      --op 1.53 \
      --ep 0.123 \
      -
else
  # The fastest but worst performance.
  # Alignments with this many seqs are always going to be
  # super gappy.
  echo "${SEQS}" \
  | mafft \
      --anysymbol \
      --inputorder \
      --quiet \
      --bl 30 \
      --parttree \
      --retree 1 \
      --maxiterate 0 \
      --nofft \
      --amino \
      --op 1.53 \
      --ep 0.123
fi
