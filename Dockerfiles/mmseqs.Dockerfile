ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

## Config variables
ENV MMSEQS_PREFIX="/opt/mmseqs"
ARG MMSEQS_VERSION="7-4e23d"
ARG MMSEQS_CMAKE_OPTIONS

ENV MPICH_PREFIX="/opt/mpich"
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

## Include MPICH
COPY --from="darcyabjones/mpich" ${MPICH_PREFIX} ${MPICH_PREFIX}

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


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV MPICH_PREFIX="/opt/mpich"
ENV MMSEQS_PREFIX="/opt/mmseqs"

COPY --from=builder ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=builder ${MMSEQS_PREFIX} ${MMSEQS_PREFIX}

RUN  apt-get update \
  && apt-get install -y \
       gawk \
       bash \
       grep \
       libstdc++6 \
       libgomp1 \
       zlib1g \
       libbz2-1.0 \
  && rm -rf -- /var/lib/apt/lists/*

ENV PATH=${PATH}:${MPICH_PREFIX}/bin:${MMSEQS_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

ENTRYPOINT ["mmseqs"]
CMD ["--help"]
