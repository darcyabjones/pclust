#!/usr/bin/env bash

SEQS="$(cat /dev/stdin | tr '\0' '\n')"
NSEQS=$(echo "${SEQS}" | grep -c "^>" )

if [ ${NSEQS} -le 1 ]
then
  echo "${SEQS}"
else
  echo "${SEQS}" | mafft --pileup --anysymbol --quiet --auto --amino --ep 0.123 -
fi
