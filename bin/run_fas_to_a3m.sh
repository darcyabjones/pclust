#!/usr/bin/env bash

TMP="tmp$$"
cat /dev/stdin | awk '!/^>/ {gsub("*", "X")} {print}' > "${TMP}_in.fasta"
reformat.pl fas a3m "${TMP}_in.fasta" "${TMP}_out.a3m" -M first > /dev/null
cat "${TMP}_out.a3m"

rm "${TMP}"{_in.fasta,_out.a3m}
