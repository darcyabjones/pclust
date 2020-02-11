ARG IMAGE="darcyabjones/base"
ARG PYTHON3_IMAGE
ARG MMSEQS_IMAGE
ARG HHSUITE_IMAGE
ARG DECIPHER_IMAGE
ARG FASTTREE_IMAGE
ARG FFDB_IMAGE

FROM "${PYTHON3_IMAGE}" as python3_builder
FROM "${MMSEQS_IMAGE}" as mmseqs_builder
FROM "${HHSUITE_IMAGE}" as hhsuite_builder
FROM "${DECIPHER_IMAGE}" as decipher_builder
FROM "${FASTTREE_IMAGE}" as fasttree_builder
FROM "${FFDB_IMAGE}" as ffdb_builder


FROM "${IMAGE}"

COPY --from=python3_builder "${APT_REQUIREMENTS_FILE}" /build/apt/python3.txt

ARG MMSEQS_TAG
ARG MMSEQS_PREFIX_ARG
ENV MMSEQS_PREFIX="${MMSEQS_PREFIX_ARG}"
LABEL mmseqs.version="${MMSEQS_TAG}"

ENV PATH="${MMSEQS_PREFIX}/bin:${PATH}"

COPY --from=mmseqs_builder "${MMSEQS_PREFIX}" "${MMSEQS_PREFIX}"
COPY --from=mmseqs_builder "${APT_REQUIREMENTS_FILE}" /build/apt/mmseqs.txt


ARG HHSUITE_TAG
ARG HHSUITE_PREFIX_ARG
ENV HHSUITE_PREFIX="${HHSUITE_PREFIX_ARG}"
ENV HHLIB="${HHSUITE_PREFIX_ARG}"
LABEL hhsuite.version="${HHSUITE_TAG}"

ENV PATH="${HHSUITE_PREFIX}/bin:${HHSUITE_PREFIX}/scripts:${PATH}"

COPY --from=hhsuite_builder "${HHSUITE_PREFIX}" "${HHSUITE_PREFIX}"
COPY --from=hhsuite_builder "${APT_REQUIREMENTS_FILE}" /build/apt/hhsuite.txt


ARG DECIPHER_VERSION
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"
LABEL decipher.version="${DECIPHER_VERSION}"

ENV R_LIBS_USER="${DECIPHER_PREFIX}:${R_LIBS_USER:-}"

COPY --from=decipher_builder "${DECIPHER_PREFIX}" "${DECIPHER_PREFIX}"
COPY --from=decipher_builder "${APT_REQUIREMENTS_FILE}" /build/apt/decipher.txt


ARG FASTTREE_VERSION
ARG FASTTREE_PREFIX_ARG
ENV FASTTREE_PREFIX="${FASTTREE_PREFIX_ARG}"
LABEL fasttree.version="${FASTTREE_VERSION}"

COPY --from=fasttree_builder "${FASTTREE_PREFIX}" "${FASTTREE_PREFIX}"
COPY --from=fasttree_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fasttree.txt

ENV PATH="${FASTTREE_PREFIX}/bin:${PATH}"


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
