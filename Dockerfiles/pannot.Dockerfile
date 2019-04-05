ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION}

RUN  if [ ! -d /usr/share/man/man1 ]; then \
       mkdir -p /usr/share/man/man1; \
     fi; \
     apt-get update \
     && apt-get install -y \
       libstdc++6 \
       libgomp1 \
       zlib1g \
       libbz2-1.0 \
       lzma \
       libsqlite3-0 \
       libffi6 \
       default-jre \
       openssl \
       gawk \
       wget \
       perl \
       python2.7 \
       python-numpy \
       python-biopython \
  && rm -rf -- /var/lib/apt/lists/* \
  && mkdir -p /tmp \
  && wget -P /usr/local/bin https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/ffdb.py \
  && wget -P /usr/local/bin https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/hh_reader.py \
  && chmod a+x /usr/local/bin/ffdb.py \
  && chmod a+x /usr/local/bin/hh_reader.py


ENV MPICH_PREFIX="/opt/mpich"
ENV PYTHON_PREFIX="/opt/python3"
ENV EMBOSS_PREFIX="/opt/EMBOSS"
ENV FFINDEX_PREFIX="/opt/ffindex"
ENV EMBOSS_PREFIX="/opt/EMBOSS"
ENV WEKA_381_PREFIX="/opt/weka-3-8-1"
ENV WEKA_3612_PREFIX="/opt/weka-3-6-12"
ENV EFFECTORP2_PREFIX="/opt/effectorp2"
ENV EFFECTORP_PREFIX="/opt/effectorp"
ENV APOPLASTP_PREFIX="/opt/apoplastp"
ENV LOCALIZER_PREFIX="/opt/localizer"

ENV SIGNALP3_PREFIX="/opt/signalp3"
ENV SIGNALP4_PREFIX="/opt/signalp4"
ENV TMHMM_PREFIX="/opt/tmhmm"
ENV CHLOROP_PREFIX="/opt/chlorop"
ENV TARGETP_PREFIX="/opt/targetp"
ENV PHOBIUS_PREFIX="/opt/phobius"

COPY --from="darcyabjones/python3" ${PYTHON_PREFIX} ${PYTHON_PREFIX}
COPY --from="darcyabjones/mpich" ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from="darcyabjones/ffindex" ${FFINDEX_PREFIX} ${FFINDEX_PREFIX}
COPY --from="darcyabjones/emboss" ${EMBOSS_PREFIX} ${EMBOSS_PREFIX}
COPY --from="darcyabjones/sperschneider" ${WEKA_381_PREFIX} ${WEKA_381_PREFIX}
COPY --from="darcyabjones/sperschneider" ${WEKA_3612_PREFIX} ${WEKA_3612_PREFIX}
COPY --from="darcyabjones/sperschneider" ${EFFECTORP2_PREFIX} ${EFFECTORP2_PREFIX}
COPY --from="darcyabjones/sperschneider" ${EFFECTORP_PREFIX} ${EFFECTORP_PREFIX}
COPY --from="darcyabjones/sperschneider" ${APOPLASTP_PREFIX} ${APOPLASTP_PREFIX}
COPY --from="darcyabjones/sperschneider" ${LOCALIZER_PREFIX} ${LOCALIZER_PREFIX}

COPY --from="darcyabjones/signalp3" ${SIGNALP3_PREFIX} ${SIGNALP3_PREFIX}
COPY --from="darcyabjones/signalp4" ${SIGNALP4_PREFIX} ${SIGNALP4_PREFIX}
COPY --from="darcyabjones/tmhmm" ${TMHMM_PREFIX} ${TMHMM_PREFIX}
COPY --from="darcyabjones/targetp" ${CHLOROP_PREFIX} ${CHLOROP_PREFIX}
COPY --from="darcyabjones/targetp" ${TARGETP_PREFIX} ${TARGETP_PREFIX}
COPY --from="darcyabjones/phobius" ${PHOBIUS_PREFIX} ${PHOBIUS_PREFIX}

ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib:${PYTHON_PREFIX}/lib:${FFINDEX_PREFIX}/lib
ENV LD_LIBRARY_PATH=${LIBRARY_PATH}
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include:${PYTHON_PREFIX}/include:${FFINDEX_PREFIX}/include

ENV PATH="${PATH}:${MPICH_PREFIX}/bin:${PYTHON_PREFIX}/bin:${FFINDEX_PREFIX}/bin"
ENV PATH="${EMBOSS_PREFIX}/bin:${EFFECTORP2_PREFIX}/Scripts:${EFFECTORP_PREFIX}/Scripts:${APOPLASTP_PREFIX}/Scripts:${LOCALIZER_PREFIX}/Scripts:${PATH}"
ENV PATH="${SIGNALP4_PREFIX}:${SIGNALP3_PREFIX}:${TMHMM_PREFIX}/bin:${CHLOROP_PREFIX}:${TARGETP_PREFIX}:${PHOBIUS_PREFIX}:${PATH}"
