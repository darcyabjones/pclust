ARG IMAGE
ARG MMSEQS_IMAGE

FROM "${MMSEQS_IMAGE}" as mmseqs_builder

FROM "${IMAGE}" as decipher_builder

ARG DECIPHER_VERSION
ARG DECIPHER_URL
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       wget \
  && add_runtime_dep \
       r-base \
       r-cran-ape \
       r-cran-phangorn \
       r-cran-rsqlite \
       r-bioc-biostrings \
       r-cran-optparse \
       r-cran-magrittr \
  && apt_install_from_file "${APT_REQUIREMENTS_FILE}" \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O "decipher.tar.gz" "${DECIPHER_URL}" \
  && tar -zxf decipher.tar.gz \
  && cd DECIPHER*/ \
  && R CMD build --no-build-vignettes --no-manual . \
  && mkdir -p "${DECIPHER_PREFIX}" \
  && R CMD INSTALL --build --library="${DECIPHER_PREFIX}" ./DECIPHER_*.tar.gz


FROM "${IMAGE}"

ARG DECIPHER_VERSION
ARG DECIPHER_PREFIX_ARG
ENV DECIPHER_PREFIX="${DECIPHER_PREFIX_ARG}"
LABEL decipher.version="${DECIPHER_VERSION}"

ENV R_LIBS_USER="${DECIPHER_PREFIX}:${R_LIBS_USER:-}"

COPY --from=decipher_builder "${DECIPHER_PREFIX}" "${DECIPHER_PREFIX}"
COPY --from=decipher_builder "${APT_REQUIREMENTS_FILE}" /build/apt/decipher.txt

ARG MMSEQS_TAG
ARG MMSEQS_PREFIX_ARG
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
