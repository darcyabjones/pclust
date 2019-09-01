ARG IMAGE

FROM "${IMAGE}" as hhsuite_builder

ARG HHSUITE_TAG
ARG HHSUITE_REPO
ARG HHSUITE_PREFIX_ARG
ARG HHSUITE_CMAKE_OPTIONS
ENV HHSUITE_PREFIX="${HHSUITE_PREFIX_ARG}"
ENV HHLIB="${HHSUITE_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       ca-certificates \
       cmake \
       git \
       libmpich-dev \
       ninja-build \
       xxd \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${HHSUITE_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${HHSUITE_TAG}" \
  && git submodule update --init \
  && cd lib/ffindex \
  && git checkout master \
  && mkdir build \
  && cd build \
  && cmake \
       -G Ninja \
       -DHAVE_MPI=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${HHSUITE_PREFIX}" \
       .. \
  && ninja \
  && ninja install \
  && cd /tmp \
  && mkdir build \
  && cd build \
  && cmake \
       -G Ninja \
       -DHAVE_MPI=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${HHSUITE_PREFIX}" \
       .. \
  && ninja \
  && ninja install \ 
  && add_runtime_dep \
       libgomp1 \
       libstdc++6 \
       mpich \
       zlib1g


FROM "${IMAGE}"

ARG HHSUITE_VERSION
ARG HHSUITE_PREFIX_ARG
ENV HHSUITE_PREFIX="${HHSUITE_PREFIX_ARG}"
ENV HHLIB="${HHSUITE_PREFIX_ARG}"
LABEL hhsuite.version="${HHSUITE_TAG}"

ENV PATH="${HHSUITE_PREFIX}/bin:${HHSUITE_PREFIX}/scripts:${PATH}"

COPY --from=hhsuite_builder "${HHSUITE_PREFIX}" "${HHSUITE_PREFIX}"
COPY --from=hhsuite_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hhsuite.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
