ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV MPICH_PREFIX="/opt/mpich"
ENV HHLIB_PREFIX="/opt/hhsuite"
ENV MMSEQS_PREFIX="/opt/mmseqs"
ENV MAFFT_PREFIX="/opt/mafft"
ENV PYTHON3_PREFIX="/opt/python3"

ENV HHLIB="${HHLIB_PREFIX}"
ENV PYTHONHOME="${PYTHON3_PREFIX}"

COPY --from=darcyabjones/mpich ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=darcyabjones/hhsuite ${HHLIB_PREFIX} ${HHLIB_PREFIX}
COPY --from=darcyabjones/mmseqs ${MMSEQS_PREFIX} ${MMSEQS_PREFIX}
COPY --from=darcyabjones/mafft ${MAFFT_PREFIX} ${MAFFT_PREFIX}
COPY --from=darcyabjones/python3 ${PYTHON3_PREFIX} ${PYTHON3_PREFIX}

RUN  apt-get update \
  && apt-get install -y \
       libbz2-1.0 \
       libffi6 \
       libgomp1 \
       libsqlite3-0 \
       libstdc++6 \
       lzma \
       openssl \
       wget \
       zlib1g \
  && rm -rf -- var/lib/apt/lists/* \
  && wget -P /usr/local/bin https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/ffdb.py \
  && wget -P /usr/local/bin https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/hh_reader.py \
  && chmod a+x /usr/local/bin/ffdb.py \
  && chmod a+x /usr/local/bin/hh_reader.py

ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib:${PYTHON3_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include:${PYTHON3_PREFIX}/include
ENV PATH=${PATH}:${MPICH_PREFIX}/bin:${HHLIB_PREFIX}/bin:${MMSEQS_PREFIX}/bin:${MAFFT_PREFIX}/bin:${PYTHON3_PREFIX}/bin
