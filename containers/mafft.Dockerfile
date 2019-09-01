ARG IMAGE
ARG MMSEQS_IMAGE

FROM "${MMSEQS_IMAGE}" as mmseqs_builder

FROM "${IMAGE}" as mafft_builder

ARG MAFFT_VERSION="7.407"
ARG MAFFT_URL="https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz"
ARG MAFFT_PREFIX_ARG="/opt/mafft/${MAFFT_VERSION}"
ENV MAFFT_PREFIX="${MAFFT_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       wget \
  && rm -rf -- /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O mafft.tar.gz "${MAFFT_URL}" \
  && tar -xf mafft.tar.gz \
  && cd mafft-*/ \
  && sed -i "s~PREFIX = /usr/local~PREFIX = ${MAFFT_PREFIX}~" core/Makefile \
  && cd core \
  && make clean \
  && make install \
  && add_runtime_dep \
       gawk \
       bash \
       grep \
       libgomp1


FROM "${IMAGE}"

ARG MAFFT_VERSION
ARG MAFFT_PREFIX_ARG="/opt/mafft/${MAFFT_VERSION}"
ENV MAFFT_PREFIX="${MAFFT_PREFIX_ARG}"
LABEL mafft.version="${MAFFT_VERSION}"

COPY --from=mafft_builder "${MAFFT_PREFIX}" "${MAFFT_PREFIX}"
COPY --from=mafft_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mafft.txt

ENV PATH="${MAFFT_PREFIX}/bin:${PATH}"


ARG MMSEQS_TAG="7-4e23d"
ARG MMSEQS_PREFIX_ARG="/opt/mmseqs/${MMSEQS_TAG}"
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt


RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
