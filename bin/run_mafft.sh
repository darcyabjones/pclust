#!/usr/bin/env bash

SEQS="$(cat /dev/stdin | tr '\0' '\n')"
NSEQS=$(echo "${SEQS}" | grep -c "^>" )

if [ ${NSEQS} -le 1 ]
then
  echo "${SEQS}"
else
  echo "${SEQS}" | mafft --retree 2 --maxiterate 2 --amino -
fi
