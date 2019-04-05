ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV FFINDEX_PREFIX="/opt/ffindex"
ARG FFINDEX_VERSION="0.9.9.9"

ENV MPICH_PREFIX="/opt/mpich"
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

COPY --from="darcyabjones/mpich" ${MPICH_PREFIX} ${MPICH_PREFIX}

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


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV FFINDEX_PREFIX="/opt/ffindex"
ENV MPICH_PREFIX="/opt/mpich"

COPY --from=builder ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=builder ${FFINDEX_PREFIX} ${FFINDEX_PREFIX}

RUN  apt-get update \
  && apt-get install -y zlib1g libbz2-1.0 \
  && rm -rf -- /var/lib/apt/lists/*

ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib:${FFINDEX_PREFIX}/lib
ENV LD_LIBRARY_PATH=${LIBRARY_PATH}
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include:${FFINDEX_PREFIX}/include

ENV PATH="${PATH}:${MPICH_PREFIX}/bin:${FFINDEX_PREFIX}/bin"
