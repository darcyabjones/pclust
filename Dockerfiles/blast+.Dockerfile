FROM frolvlad/alpine-glibc:alpine-3.8_glibc-2.28
# Unfortunately blast relies on some backtracking features that musl doesn't have
# Need to use glibc

LABEL maintainer="darcy.ab.jones@gmail.com"
LABEL version="0.1"

RUN apk add --no-cache libbz2 libidn libgomp bash

# Setup blast
ENV BLAST_PREFIX=/opt/blast
ENV BLAST_VERSION=2.7.1
ENV NCBI_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+

WORKDIR /opt
RUN    wget ${NCBI_URL}/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && wget ${NCBI_URL}/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz.md5 \
    && md5sum -c ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz.md5 \
    && tar -zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz* \
    && mv ncbi-blast-${BLAST_VERSION}+ ${BLAST_PREFIX}

ENV BLASTDB=/blastdb
ENV PATH=${BLAST_PREFIX}/bin:${PATH}
