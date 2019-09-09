ARG IMAGE
ARG FFDB_IMAGE

FROM "${FFDB_IMAGE}" as ffdb_builder

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
  && mkdir build \
  && cd build \
  && cmake \
       -G Ninja \
       -DHAVE_MPI=1 \
       -DHAVE_SSE2=1 \
       -DHAVE_AVX2=0 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${HHSUITE_PREFIX}" \
       .. \
  && ninja \
  && ninja install \
  && mv "${HHSUITE_PREFIX}/bin" "${HHSUITE_PREFIX}/bin-sse" \
  && cd /tmp \
  && rm -rf -- build \
  && mkdir build \
  && cd build \
  && cmake \
       -G Ninja \
       -DHAVE_MPI=1 \
       -DHAVE_AVX2=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${HHSUITE_PREFIX}" \
       .. \
  && ninja \
  && ninja install \
  && mv "${HHSUITE_PREFIX}/bin" "${HHSUITE_PREFIX}/bin-avx2" \
  && mkdir "${HHSUITE_PREFIX}/bin" \
  && echo '#!/usr/bin/env bash'                                  > "${HHSUITE_PREFIX}/template" \
  && echo 'EXE=$(basename $0)'                                  >> "${HHSUITE_PREFIX}/template" \
  && echo 'if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then' >> "${HHSUITE_PREFIX}/template" \
  && echo '    exec "${HHSUITE_PREFIX}/bin-avx2/${EXE}" "$@"'   >> "${HHSUITE_PREFIX}/template" \
  && echo 'else'                                                >> "${HHSUITE_PREFIX}/template" \
  && echo '    exec "${HHSUITE_PREFIX}/bin-sse/${EXE}" "$@"'    >> "${HHSUITE_PREFIX}/template" \
  && echo 'fi'                                                  >> "${HHSUITE_PREFIX}/template" \
  && chmod a+x "${HHSUITE_PREFIX}/template" \
  && for f in ${HHSUITE_PREFIX}/bin-sse/*; do \
       cp "${HHSUITE_PREFIX}/template" ${HHSUITE_PREFIX}/bin/$(basename "${f}"); \
     done \
  && rm "${HHSUITE_PREFIX}/template" \
  && add_runtime_dep \
       libgomp1 \
       libstdc++6 \
       mpich \
       zlib1g


FROM "${IMAGE}"

ARG HHSUITE_TAG
ARG HHSUITE_PREFIX_ARG
ENV HHSUITE_PREFIX="${HHSUITE_PREFIX_ARG}"

ENV HHLIB="${HHSUITE_PREFIX_ARG}"
LABEL hhsuite.version="${HHSUITE_TAG}"


ENV PATH="${HHSUITE_PREFIX}/bin:${HHSUITE_PREFIX}/scripts:${PATH}"

COPY --from=hhsuite_builder "${HHSUITE_PREFIX}" "${HHSUITE_PREFIX}"
COPY --from=hhsuite_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hhsuite.txt


ARG FFDB_TAG
ARG FFDB_PREFIX_ARG="/opt/ffdb/${FFDB_TAG}"
ENV FFDB_PREFIX="${FFDB_PREFIX_ARG}"
LABEL ffdb.version="${FFDB_TAG}"

ENV PATH "${FFDB_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${FFDB_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=ffdb_builder "${FFDB_PREFIX}" "${FFDB_PREFIX}"
COPY --from=ffdb_builder "${APT_REQUIREMENTS_FILE}" /build/apt/ffdb.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
