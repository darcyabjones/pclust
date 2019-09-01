#!/usr/bin/env bash

TMP="tmp$$"
cat /dev/stdin > "${TMP}_in.fasta"

reformat.pl fas a3m "${TMP}_in.fasta" "${TMP}_out.a3m" -M 100 > /dev/null
cat "${TMP}_out.a3m"

rm "${TMP}"{_in.fasta,_out.a3m}
