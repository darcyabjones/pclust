FROM alpine:3.8 as psibuilder

RUN apk add --no-cache build-base

ENV PSIPRED_VERSION=4.02
ENV PSIPRED_PREFIX=/opt/psipred
ENV PSIPRED_URL=http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred

WORKDIR /opt
RUN    wget ${PSIPRED_URL}/psipred.${PSIPRED_VERSION}.tar.gz
RUN tar -zxf psipred.${PSIPRED_VERSION}.tar.gz \
    && rm psipred.${PSIPRED_VERSION}.tar.gz \
    && cd ${PSIPRED_PREFIX}/src \
    && make \
    && make install


FROM frolvlad/alpine-glibc:alpine-3.8_glibc-2.28


RUN apk add --no-cache \
      bash \
      libbz2 \
      libidn \
      tcsh \
      grep \
      libstdc++ \
      libgomp \
      python \
      perl

# Include psipred
COPY --from=psibuilder /opt/psipred /opt/psipred
ENV PSIPRED_PREFIX=/opt/psipred

# Include MPICH
COPY --from=pclust/hhblits:latest /opt/mpich /opt/mpich
ENV MPICH_PREFIX=/opt/mpich
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

# Include HHblits
COPY --from=pclust/hhblits:latest /opt/hh-suite /opt/hh-suite

ENV HHLIB=/opt/hh-suite
ENV PATH="/opt/hh-suite/bin:/opt/hh-suite/scripts:${PATH}"

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

# Modify psipred path
RUN    sed -i "s~ncbidir = /usr/local/bin~ncbidir = ${BLAST_PREFIX}/bin~g" ${PSIPRED_PREFIX}/BLAST+/runpsipredplus \
    && sed -i "s~/cluster/toolkit/production/bioprogs~/opt~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/pdb/all~/data/pdb~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/pdb/all~/data/pdb~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/pdb/all~/data/pdb~" ${HHLIB}/scripts/HHPaths.pm \


ENV BLASTDB=/blast
ENV PATH=${BLAST_PREFIX}/bin:${PSIPRED_PREFIX}/bin:${PSIPRED_PREFIX}:${PSIPRED_PREFIX}/BLAST+:${PATH}
