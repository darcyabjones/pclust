from alpine:3.8 as mpibuilder
RUN apk add --no-cache build-base

# Source is available at http://www.mpich.org/static/downloads/

# Build Options:
# See installation guide of target MPICH version
# Ex: http://www.mpich.org/static/downloads/3.2/mpich-3.2-installguide.pdf

ENV MPICH_PREFIX="/opt/mpich"
ARG MPICH_VERSION="3.2.1"
ARG MPICH_CONFIGURE_OPTIONS="--disable-fortran --prefix=${MPICH_PREFIX}"
ARG MPICH_MAKE_OPTIONS

# Download, build, and install MPICH
WORKDIR /tmp/mpich-src
RUN wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
      && tar xfz mpich-${MPICH_VERSION}.tar.gz  \
      && cd mpich-${MPICH_VERSION}  \
      && ./configure ${MPICH_CONFIGURE_OPTIONS}  \
      && make ${MPICH_MAKE_OPTIONS} \
      && make install

RUN rm -rf /tmp/*

from alpine:3.8 as hhbuilder

# Include MPICH
COPY --from=mpibuilder /opt/mpich /opt/mpich
ENV MPICH_PREFIX=/opt/mpich
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

RUN apk add --no-cache build-base cmake musl-dev git ninja

WORKDIR /opt/hh-suite
RUN  git clone https://github.com/soedinglab/hh-suite.git . \
  && git fetch --tags \
  && git checkout tags/v3.0-beta.3 \
  && git submodule update --init

WORKDIR /opt/hh-suite/lib/ffindex
RUN git checkout master

WORKDIR /opt/hh-suite/build
RUN cmake \
      -G Ninja \
      -DHAVE_SSE2=1 \
      -DHAVE_MPI=1 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/opt/hh-suite ..

RUN ninja && ninja install

WORKDIR /opt/hh-suite
RUN rm -rf /opt/hh-suite/build

FROM soedinglab/mmseqs2:version-5

MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache bash grep libstdc++ libgomp

# Include MPICH
COPY --from=mpibuilder /opt/mpich /opt/mpich
ENV MPICH_PREFIX=/opt/mpich
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

# Include HHblits
COPY --from=hhbuilder /opt/hh-suite /opt/hh-suite

ENV HHLIB=/opt/hh-suite
ENV PATH="/opt/hh-suite/bin:${PATH}"
