FROM frolvlad/alpine-glibc:alpine-3.8_glibc-2.28

LABEL maintainer="darcy.ab.jones@gmail.com"
LABEL version="0.1"

RUN apk add --no-cache bash wget

# Setup blast
ENV BLAST_PREFIX=/opt/blast
ENV BASE_URL=https://anaconda.org/bioconda/blast-legacy/2.2.26/download/linux-64/blast-legacy-2.2.26-1.tar.bz2

WORKDIR ${BLAST_PREFIX}
RUN    wget ${BASE_URL} \
    && tar -xjf blast-legacy-2.2.26-1.tar.bz2 \
    && rm blast-legacy-2.2.26-1.tar.bz2 \
    && rm -rf -- info \
    && wget -r ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/ \
    && mv ftp.ncbi.nlm.nih.gov/blast/matrices ${BLAST_PREFIX}/data \
    && rm -rf -- ftp.ncbi.nlm.nih.gov \
    && echo -e "[NCBI]\nDATA=${BLAST_PREFIX}/data" > /root/.ncbirc


ENV BLASTDB=/blastdb
ENV PATH=${BLAST_PREFIX}/bin:${PATH}
WORKDIR /
