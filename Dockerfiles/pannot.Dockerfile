ARG DEBIAN_VERSION="stretch-20190228-slim"

# INSTALL MPICH
# =============

FROM debian:${DEBIAN_VERSION} as mpibuilder
RUN  apt-get update \
  && apt-get install -y build-essential wget\
  && rm -rf -- /var/lib/apt/lists/*

# Source is available at http://www.mpich.org/static/downloads/

# Build Options:
# See installation guide of target MPICH version
# Ex: http://www.mpich.org/static/downloads/3.2/mpich-3.2-installguide.pdf

## Config variables
ENV MPICH_PREFIX="/opt/mpich"
ARG MPICH_VERSION="3.2.1"
ARG MPICH_CONFIGURE_OPTIONS="--disable-fortran --prefix=${MPICH_PREFIX}"
ARG MPICH_MAKE_OPTIONS

## Download, build, and install MPICH
WORKDIR /tmp/mpich-src
RUN  wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
  && tar xfz mpich-${MPICH_VERSION}.tar.gz  \
  && cd mpich-${MPICH_VERSION}  \
  && ./configure ${MPICH_CONFIGURE_OPTIONS}  \
  && make ${MPICH_MAKE_OPTIONS} \
  && make install \
  && rm -rf /tmp/*


# INSTALL PYTHON
# ==============

# Debians packaged python 3 is still 3.5, which is too old.

FROM debian:${DEBIAN_VERSION} as pybuilder

## Config variables
ENV PYTHON_PREFIX="/opt/python"
ARG PYTHON3_VERSION="3.7.2"

RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
       liblzma-dev \
       zlib1g-dev \
       libbz2-dev \
       libsqlite3-0 \
       libsqlite3-dev \
       libffi-dev \
       libssl1.0 \
       libssl1.0-dev \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/python3
RUN  wget https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tgz \
  && tar xfz Python-${PYTHON3_VERSION}.tgz \
  && cd Python-${PYTHON3_VERSION} \
  && ./configure --enable-optimizations --prefix ${PYTHON_PREFIX} \
  && make \
  && make install


# INSTALL EMBOSS
# ==============

FROM debian:${DEBIAN_VERSION} as embossbuilder

## Config variables
ENV EMBOSS_PREFIX="/opt/EMBOSS"
ARG EMBOSS_VERSION="6.5.7"
ARG EMBOSS_URL="ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz"

RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/emboss
RUN  wget ${EMBOSS_URL} \
  && tar xfz EMBOSS-${EMBOSS_VERSION}.tar.gz \
  && cd EMBOSS-${EMBOSS_VERSION} \
  && ./configure --without-x --prefix ${EMBOSS_PREFIX} \
  && make \
  && make install


# INSTALL FFINDEX
# ===============

FROM debian:${DEBIAN_VERSION} as ffbuilder

ENV FFINDEX_PREFIX="/opt/ffindex"
ARG FFINDEX_VERSION="0.9.9.9"

ENV MPICH_PREFIX="/opt/mpich"
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

## Include MPICH
COPY --from=mpibuilder ${MPICH_PREFIX} ${MPICH_PREFIX}

RUN  apt-get update \
  && apt-get install -y build-essential git zlib1g-dev libbz2-dev \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/ffindex
RUN  git clone https://github.com/ahcm/ffindex.git . \
  && git fetch --tags \
  && git checkout tags/${FFINDEX_VERSION} \
  && make HAVE_MPI=1 \
  && make test \
  && make install INSTALL_DIR=${FFINDEX_PREFIX}


# INSTALL WEKA
# ============

# Weka is installed along with effectorp etc.


# INSTALL SPERSCHNEIDER TOOLS
# ===========================

FROM debian:${DEBIAN_VERSION} as sperbuilder

ENV EMBOSS_PREFIX="/opt/EMBOSS"
ENV WEKA_381_PREFIX="/opt/weka-3-8-1"
ENV WEKA_3612_PREFIX="/opt/weka-3-6-12"
ENV EFFECTORP2_PREFIX="/opt/effectorp2"
ENV EFFECTORP_PREFIX="/opt/effectorp"
ENV APOPLASTP_PREFIX="/opt/apoplastp"
ENV LOCALIZER_PREFIX="/opt/localizer"

RUN  apt-get update \
  && apt-get install -y wget unzip \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /opt
RUN  wget http://effectorp.csiro.au/EffectorP_1.0.tar.gz \
  && tar zxf EffectorP_1.0.tar.gz \
  && rm EffectorP_1.0.tar.gz \
  && mv EffectorP_1.0 ${EFFECTORP_PREFIX} \
  && cd ${EFFECTORP_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && unzip -d /opt weka-3-6-12.zip \
  && rm weka-3-6-12.zip \
  && chmod a+x EffectorP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" EffectorP.py \
  && ln -sf "${WEKA_3612_PREFIX}" "${EFFECTORP_PREFIX}/Scripts/weka-3-6-12" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" EffectorP.py \
  && ln -sf ${PWD}/EffectorP.py ${PWD}/EffectorP-1.0.py

WORKDIR /opt
RUN  wget http://localizer.csiro.au/LOCALIZER_1.0.4.tar.gz \
  && tar zxf LOCALIZER_1.0.4.tar.gz \
  && rm LOCALIZER_1.0.4.tar.gz \
  && mv LOCALIZER_1.0.4 ${LOCALIZER_PREFIX} \
  && cd ${LOCALIZER_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && rm weka-3-6-12.zip \
  && chmod a+x LOCALIZER.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" LOCALIZER.py \
  && ln -sf "${WEKA_3612_PREFIX}" "${LOCALIZER_PREFIX}/Scripts/weka-3-6-12" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" LOCALIZER.py

WORKDIR /opt
RUN  wget http://apoplastp.csiro.au/ApoplastP_1.0.1.tar.gz \
  && tar zxf ApoplastP_1.0.1.tar.gz \
  && rm ApoplastP_1.0.1.tar.gz \
  && mv ApoplastP_1.0.1 ${APOPLASTP_PREFIX} \
  && cd ${APOPLASTP_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && unzip -d /opt weka-3-8-1.zip \
  && rm weka-3-8-1.zip \
  && chmod a+x ApoplastP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" ApoplastP.py \
  && ln -sf "${WEKA_381_PREFIX}" "${APOPLASTP_PREFIX}/Scripts/weka-3-8-1" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" ApoplastP.py

WORKDIR /opt
RUN  wget http://effectorp.csiro.au/EffectorP_2.0.tar.gz \
  && tar zxf EffectorP_2.0.tar.gz \
  && rm EffectorP_2.0.tar.gz \
  && mv EffectorP_2.0 ${EFFECTORP2_PREFIX} \
  && cd ${EFFECTORP2_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && rm weka-3-8-1.zip \
  && chmod a+x EffectorP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" EffectorP.py \
  && ln -sf "${WEKA_381_PREFIX}" "${EFFECTORP2_PREFIX}/Scripts/weka-3-8-1" \
  && mkdir "${EFFECTORP2_PREFIX}/Scripts/EMBOSS-6.5.7" \
  && ln -sf "${EMBOSS_PREFIX}/bin" "${EFFECTORP2_PREFIX}/Scripts/EMBOSS-6.5.7/emboss" \
  && ln -sf "${PWD}/EffectorP.py" "${PWD}/EffectorP-2.0.py"


# INSTALL CBS TOOLS
# =================

FROM debian:${DEBIAN_VERSION} as cbsbuilder

ENV SIGNALP3_PREFIX="/opt/signalp3"
ENV SIGNALP4_PREFIX="/opt/signalp4"
ENV TMHMM_PREFIX="/opt/tmhmm"
ENV CHLOROP_PREFIX="/opt/chlorop"
ENV TARGETP_PREFIX="/opt/targetp"
ENV PHOBIUS_PREFIX="/opt/phobius"

COPY signalp-3.0.Linux.tar.Z /opt/signalp-3.0.Linux.tar.Z
COPY signalp-4.1f.Linux.tar.gz /opt/signalp-4.1f.Linux.tar.gz
COPY tmhmm-2.0c.Linux.tar.gz /opt/tmhmm-2.0c.Linux.tar.gz
COPY chlorop-1.1.Linux.tar.Z /opt/chlorop-1.1.Linux.tar.Z
COPY targetp-1.1b.Linux.tar.Z /opt/targetp-1.1b.Linux.tar.Z
COPY phobius101_linux.tar.gz /opt/phobius101_linux.tar.gz


WORKDIR /opt
RUN  tar xf signalp-3.0.Linux.tar.Z \
  && rm signalp-3.0.Linux.tar.Z \
  && mv signalp-3.0 ${SIGNALP3_PREFIX} \
  && sed -i s~SIGNALP=/usr/opt/signalp-3.0~SIGNALP=${SIGNALP3_PREFIX}~ ${SIGNALP3_PREFIX}/signalp \
  && sed -i s~AWK=/usr/bin/gawk~AWK="/usr/bin/env gawk"~ ${SIGNALP3_PREFIX}/signalp \
  && cp ${SIGNALP3_PREFIX}/signalp ${SIGNALP3_PREFIX}/signalp-3.0


WORKDIR /opt
RUN  tar xf signalp-4.1f.Linux.tar.gz \
  && rm signalp-4.1f.Linux.tar.gz \
  && mv signalp-4.1 ${SIGNALP4_PREFIX} \
  && sed -i s~/usr/cbs/bio/src/signalp-4.1~${SIGNALP4_PREFIX}~ ${SIGNALP4_PREFIX}/signalp \
  && sed -i s~/var/tmp~/tmp~ ${SIGNALP4_PREFIX}/signalp \
  && sed -i s~MAX_ALLOWED_ENTRIES=10000~MAX_ALLOWED_ENTRIES=999999999~ ${SIGNALP4_PREFIX}/signalp \
  && cp ${SIGNALP4_PREFIX}/signalp ${SIGNALP4_PREFIX}/signalp-4.1


WORKDIR /opt
RUN  tar xf tmhmm-2.0c.Linux.tar.gz \
  && rm tmhmm-2.0c.Linux.tar.gz \
  && mv tmhmm-2.0c ${TMHMM_PREFIX} \
  && sed -i s~/usr/local/bin/perl~"/usr/bin/env perl"~ ${TMHMM_PREFIX}/bin/tmhmm \
  && sed -i s~/usr/local/bin/perl~"/usr/bin/env perl"~ ${TMHMM_PREFIX}/bin/tmhmmformat.pl


WORKDIR /opt
RUN  tar xf chlorop-1.1.Linux.tar.Z \
  && rm chlorop-1.1.Linux.tar.Z \
  && mv chlorop-1.1 ${CHLOROP_PREFIX} \
  && sed -i s~/usr/cbs/packages/chlorop/currdist/chlorop-1.1~${CHLOROP_PREFIX}~ ${CHLOROP_PREFIX}/chlorop


WORKDIR /opt
RUN  tar xf targetp-1.1b.Linux.tar.Z \
  && rm targetp-1.1b.Linux.tar.Z \
  && mv targetp-1.1 ${TARGETP_PREFIX} \
  && sed -i s~/usr/cbs/packages/targetp/currdist/targetp-1.1~${TARGETP_PREFIX}~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/scratch~/tmp~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/bin/perl~"/usr/bin/env perl"~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~AWK=/usr/freeware/bin/gawk~AWK="/usr/bin/env gawk"~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/chlorop~${CHLOROP_PREFIX}/chlorop~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/signalp~${SIGNALP3_PREFIX}/signalp~ ${TARGETP_PREFIX}/targetp


WORKDIR /opt
RUN  tar xf phobius101_linux.tar.gz \
  && rm phobius101_linux.tar.gz \
  && mv tmp/tmpTB6e6P/phobius ${PHOBIUS_PREFIX} \
  && sed -i "s~PHOBIUS_DIR/decodeanhmm~PHOBIUS_DIR/decodeanhmm.64bit~g" ${PHOBIUS_PREFIX}/phobius.pl \
  && sed -i '181s/predstr/predstr=""/' ${PHOBIUS_PREFIX}/phobius.pl \
  && sed -i '244a \$predstr="";' ${PHOBIUS_PREFIX}/phobius.pl \
  && ln -sf ${PHOBIUS_PREFIX}/phobius.pl ${PHOBIUS_PREFIX}/phobius


# FINAL IMAGE
# ===========
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
ENV PYTHON_PREFIX="/opt/python"
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

COPY --from=mpibuilder ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=pybuilder ${PYTHON_PREFIX} ${PYTHON_PREFIX}
COPY --from=embossbuilder ${EMBOSS_PREFIX} ${EMBOSS_PREFIX}
COPY --from=ffbuilder ${FFINDEX_PREFIX} ${FFINDEX_PREFIX}
COPY --from=sperbuilder ${WEKA_381_PREFIX} ${WEKA_381_PREFIX}
COPY --from=sperbuilder ${WEKA_3612_PREFIX} ${WEKA_3612_PREFIX}
COPY --from=sperbuilder ${EFFECTORP2_PREFIX} ${EFFECTORP2_PREFIX}
COPY --from=sperbuilder ${EFFECTORP_PREFIX} ${EFFECTORP_PREFIX}
COPY --from=sperbuilder ${APOPLASTP_PREFIX} ${APOPLASTP_PREFIX}
COPY --from=sperbuilder ${LOCALIZER_PREFIX} ${LOCALIZER_PREFIX}

COPY --from=cbsbuilder ${SIGNALP3_PREFIX} ${SIGNALP3_PREFIX}
COPY --from=cbsbuilder ${SIGNALP4_PREFIX} ${SIGNALP4_PREFIX}
COPY --from=cbsbuilder ${TMHMM_PREFIX} ${TMHMM_PREFIX}
COPY --from=cbsbuilder ${CHLOROP_PREFIX} ${CHLOROP_PREFIX}
COPY --from=cbsbuilder ${TARGETP_PREFIX} ${TARGETP_PREFIX}
COPY --from=cbsbuilder ${PHOBIUS_PREFIX} ${PHOBIUS_PREFIX}

ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib:${PYTHON_PREFIX}/lib:${FFINDEX_PREFIX}/lib
ENV LD_LIBRARY_PATH=${LIBRARY_PATH}
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include:${PYTHON_PREFIX}/include:${FFINDEX_PREFIX}/include

ENV PATH="${PATH}:${MPICH_PREFIX}/bin:${PYTHON_PREFIX}/bin:${FFINDEX_PREFIX}/bin"
ENV PATH="${EMBOSS_PREFIX}/bin:${EFFECTORP2_PREFIX}/Scripts:${EFFECTORP_PREFIX}/Scripts:${APOPLASTP_PREFIX}/Scripts:${LOCALIZER_PREFIX}/Scripts:${PATH}"
ENV PATH="${SIGNALP4_PREFIX}:${SIGNALP3_PREFIX}:${TMHMM_PREFIX}/bin:${CHLOROP_PREFIX}:${TARGETP_PREFIX}:${PHOBIUS_PREFIX}:${PATH}"

ENV PYTHONHOME="/usr:/usr/local:${PYTHON_PREFIX}:${PYTHONHOME}"
ENV PYTHONPATH="/usr/lib/python2.7/lib-dynload/:/usr/lib/python2.7/:/usr/bin:/usr/local/bin:${PYTHON_PREFIX}/bin:${PYTHON_PREFIX}/lib:${PYTHONPATH}"
