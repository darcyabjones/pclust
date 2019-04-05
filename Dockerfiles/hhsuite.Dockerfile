ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV MPICH_PREFIX=/opt/mpich
ENV PATH=${PATH}:${MPICH_PREFIX}/bin
ENV LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib
ENV CPATH=${CPATH}:${MPICH_PREFIX}/include

# Include MPICH
COPY --from=darcyabjones/mpich ${MPICH_PREFIX} ${MPICH_PREFIX}

RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       cmake \
       git \
       ninja-build \
       xxd \
  && rm -rf /var/lib/apt/lists/*

ENV HHSUITE_PREFIX="/opt/hhsuite"
ARG HHSUITE_VERSION="v3.1.0"
ARG HHSUITE_CMAKE_OPTIONS

WORKDIR /tmp
RUN  git clone https://github.com/soedinglab/hh-suite.git . \
  && git fetch --tags \
  && git checkout tags/v3.2.0 \
  && git submodule update --init \
  && cd lib/ffindex \
  && git checkout master

WORKDIR /tmp/build
RUN cmake \
      -G Ninja \
      -DHAVE_MPI=1 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=${HHSUITE_PREFIX} .. \
  && ninja \
  && ninja install


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y \
       libstdc++6 \
       libgomp1 \
  && rm -rf /var/lib/apt/lists/*

ENV HHSUITE_PREFIX="/opt/hhsuite"
ENV MPICH_PREFIX="/opt/mpich"

COPY --from=darcyabjones/mpich ${MPICH_PREFIX} ${MPICH_PREFIX}
COPY --from=builder ${HHSUITE_PREFIX} ${HHSUITE_PREFIX}

ENV HHLIB="${HHSUITE_PREFIX}"
ENV PATH="${HHSUITE_PREFIX}/bin:${HHSUITE_PREFIX}/scripts:${MPICH_PREFIX}/bin:${PATH}"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${MPICH_PREFIX}/lib"
ENV CPATH="${CPATH}:${MPICH_PREFIX}/include"

ENTRYPOINT ["hhblits"]
CMD ["-help"]
