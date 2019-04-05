ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder
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

## Clean up the mpich stuff.
FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV MPICH_PREFIX="/opt/mpich"

COPY --from=builder ${MPICH_PREFIX} ${MPICH_PREFIX}

ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

ENTRYPOINT ["mpirun"]
CMD ["--help"]
