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


FROM pclust/blast:latest

RUN apk add --no-cache \
      bash \
      tcsh \
      grep \
      libstdc++ \
      libgomp \
      python \
      perl \
      libbz2 \
      boost-iostreams \
      boost-program_options \
      boost-thread

# Include psipred
COPY --from=psibuilder /opt/psipred /opt/psipred
ENV PSIPRED_PREFIX=/opt/psipred
ENV PATH=${PSIPRED_PREFIX}:${PATH}  

# Include xssp
COPY --from=pclust/xssp:latest /opt/xssp /opt/xssp
ENV PATH=/opt/xssp/bin:${PATH}  

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
ENV BLASTDB=/data/blast
ENV BLAST_PREFIX=/opt/blast

# Modify psipred path
RUN    sed -i "s~/cluster/toolkit/production/bioprogs~/opt~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/pdb/all~/data/pdb~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/dssp/data~/data/dssp~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~/cluster/databases/dssp/bin/dsspcmbi~/opt/xssp/bin~" ${HHLIB}/scripts/HHPaths.pm \
    && sed -i "s~ncbidir = /usr/local/bin~ncbidir = ${BLAST_PREFIX}/bin~" ${PSIPRED_PREFIX}/runpsipred
