ARG DEBIAN_VERSION="stretch-20190228-slim"

# INSTALL MPICH
# =============
FROM debian:${DEBIAN_VERSION} as mpibuilder
RUN  apt-get update \
  && apt-get install -y build-essential wget \
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


# INSTALL HHSUITE
# ===============
FROM debian:${DEBIAN_VERSION} as hhbuilder

## Config variables
ENV HHLIB_PREFIX="/opt/hh-suite"
ARG HHLIB_VERSION="v3.1.0"
ARG HHLIB_CMAKE_OPTIONS

ENV MPICH_PREFIX="/opt/mpich"
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

## Include MPICH
COPY --from=mpibuilder ${MPICH_PREFIX} ${MPICH_PREFIX}

RUN  apt-get update \
  && apt-get install -y build-essential cmake xxd ninja-build git \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/hh-suite
RUN  git clone https://github.com/soedinglab/hh-suite.git . \
  && git fetch --tags \
  && git checkout tags/${HHLIB_VERSION} \
  && git submodule update --init

WORKDIR /tmp/hh-suite/lib/ffindex
RUN git checkout master

WORKDIR /tmp/hh-suite/build
RUN  cmake \
       -G Ninja \
       ${HHLIB_CMAKE_OPTIONS} \
       -DHAVE_MPI=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX=${HHLIB_PREFIX} .. \
  && ninja \
  && ninja install


# INSTALL MAFFT
# =============
FROM debian:${DEBIAN_VERSION} as mafftbuilder

## Config variables
ENV MAFFT_PREFIX="/opt/mafft"
ARG MAFFT_VERSION="7.407"

RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf -- /var/lib/apt/lists/*

# Download, build, and install MPICH
WORKDIR /tmp/mafft
RUN  wget https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
  && tar zxf mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
  && cd mafft-${MAFFT_VERSION}-without-extensions \
  && sed -i "s~PREFIX = /usr/local~PREFIX = ${MAFFT_PREFIX}~" core/Makefile \
  && cd core \
  && make clean \
  && make install


# INSTALL MMSEQS2
# ===============
FROM debian:${DEBIAN_VERSION} as mmbuilder

## Config variables
ENV MMSEQS_PREFIX="/opt/mmseqs"
ARG MMSEQS_VERSION="7-4e23d"
ARG MMSEQS_CMAKE_OPTIONS

ENV MPICH_PREFIX="/opt/mpich"
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

## Include MPICH
COPY --from=mpibuilder ${MPICH_PREFIX} ${MPICH_PREFIX}

RUN  apt-get update \
  && apt-get install -y build-essential cmake xxd git zlib1g-dev libbz2-dev \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/mmseqs
RUN  git clone https://github.com/soedinglab/MMseqs2.git . \
  && git fetch --tags \
  && git checkout tags/${MMSEQS_VERSION} \
  && git submodule update --init

WORKDIR /tmp/mmseqs/build
RUN  cmake \
       ${MMSEQS_CMAKE_OPTIONS} \
       -DHAVE_MPI=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX=${MMSEQS_PREFIX} .. \
  && make \
  && make install


# INSTALL PYTHON3
# ===============
FROM debian:${DEBIAN_VERSION} as pybuilder

## Config variables
ENV PYTHON3_PREFIX="/opt/python3"
ARG PYTHON3_VERSION="3.7.2"

RUN  apt-get update \
  && apt-get install -y build-essential wget liblzma-dev zlib1g-dev libbz2-dev libsqlite3-0 libsqlite3-dev libffi-dev libssl1.0 libssl1.0-dev \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/python3
RUN  wget https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tgz \
  && tar xfz Python-${PYTHON3_VERSION}.tgz \
  && cd Python-${PYTHON3_VERSION} \
  && ./configure --enable-optimizations --prefix ${PYTHON3_PREFIX} \
  && make \
  && make install

# && make test \

# FINAL IMAGE
# ===========
FROM debian:${DEBIAN_VERSION}

ENV MPICH_PREFIX="/opt/mpich"
ENV HHLIB_PREFIX="/opt/hh-suite"
ENV MMSEQS_PREFIX="/opt/mmseqs"
ENV MAFFT_PREFIX="/opt/mafft"
ENV PYTHON3_PREFIX="/opt/python3"

ENV HHLIB="${HHLIB_PREFIX}"
ENV PYTHONHOME="${PYTHON3_PREFIX}"

COPY --from=mpibuilder ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=hhbuilder ${HHLIB_PREFIX} ${HHLIB_PREFIX}
COPY --from=mmbuilder ${MMSEQS_PREFIX} ${MMSEQS_PREFIX}
COPY --from=mafftbuilder ${MAFFT_PREFIX} ${MAFFT_PREFIX}
COPY --from=pybuilder ${PYTHON3_PREFIX} ${PYTHON3_PREFIX}

RUN  apt-get update \
  && apt-get install -y wget libstdc++6 libgomp1 zlib1g libbz2-1.0 lzma libsqlite3-0 libffi6 openssl \
  && rm -rf -- var/lib/apt/lists/* \
  && wget -P /usr/local/bin https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/ffdb.py \
  && wget -P https://raw.githubusercontent.com/darcyabjones/pclust/master/bin/hh_reader.py \
  && chmod a+x /usr/local/bin/ffdb.py \
  && chmod a+x /usr/local/bin/hh_reader.py

ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib:${PYTHON3_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include:${PYTHON3_PREFIX}/include
ENV PATH=${PATH}:${MPICH_PREFIX}/bin:${HHLIB_PREFIX}/bin:${MMSEQS_PREFIX}/bin:${MAFFT_PREFIX}/bin:${PYTHON3_PREFIX}/bin

ENV PYTHONHOME="${PYTHON3_PREFIX}"
ENV PYTHONPATH="${PYTHON3_PREFIX}/bin"
