ARG IMAGE
ARG FFDB_IMAGE

FROM "${FFDB_IMAGE}" as ffdb_builder

FROM "${IMAGE}" as mmseqs_builder

## Config variables
ARG MMSEQS_TAG="9-d36de"
ARG MMSEQS_CMAKE_OPTIONS=""
ARG MMSEQS_REPO="https://github.com/soedinglab/MMseqs2.git"
ARG MMSEQS_PREFIX_ARG="/opt/mmseqs/${MMSEQS_TAG}"
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       git \
       libbz2-dev \
       libmpich-dev \
       xxd \
       zlib1g-dev \
  && rm -rf -- /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${MMSEQS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${MMSEQS_TAG}" \
  && git submodule update --init \
  && mkdir -p build \
  && cd build \
  && cmake \
       ${MMSEQS_CMAKE_OPTIONS} \
       -DHAVE_MPI=1 \
       -DHAVE_AVX2=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${MMSEQS_PREFIX}" .. \
  && make \
  && make install \
  && mv "${MMSEQS_PREFIX}/bin/mmseqs" "${MMSEQS_PREFIX}/bin/mmseqs.avx2" \
  && cmake \
       ${MMSEQS_CMAKE_OPTIONS} \
       -DHAVE_MPI=1 \
       -DHAVE_AVX2=0 \
       -DHAVE_SSE4_1=1 \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_INSTALL_PREFIX="${MMSEQS_PREFIX}" .. \
  && make \
  && make install \
  && mv "${MMSEQS_PREFIX}/bin/mmseqs" "${MMSEQS_PREFIX}/bin/mmseqs.sse4" \
  && cp ../util/mmseqs_wrapper.sh "${MMSEQS_PREFIX}/bin/mmseqs" \
  && sed -i 's~/usr/local/bin/mmseqs_avx2~"${MMSEQS_PREFIX}/bin/mmseqs.avx2"~g' "${MMSEQS_PREFIX}/bin/mmseqs" \
  && sed -i 's~/usr/local/bin/mmseqs_sse42~"${MMSEQS_PREFIX}/bin/mmseqs.sse4"~g' "${MMSEQS_PREFIX}/bin/mmseqs" \
  && add_runtime_dep \
       gawk \
       bash \
       grep \
       libbz2-1.0 \
       libgomp1 \
       libstdc++6 \
       mpich \
       zlib1g


FROM "${IMAGE}"

ARG MMSEQS_TAG="7-4e23d"
ARG MMSEQS_PREFIX_ARG="/opt/mmseqs/${MMSEQS_TAG}"
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt


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
