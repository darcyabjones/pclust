ARG IMAGE
ARG HHSUITE_IMAGE

FROM "${HHSUITE_IMAGE}" as hhsuite_builder

FROM "${IMAGE}" as ffdb_builder

ARG FFDB_TAG
ARG FFDB_REPO="https://github.com/darcyabjones/ffdb.git"
ARG FFDB_PREFIX_ARG="/opt/ffdb/${FFDB_TAG}"
ENV FFDB_PREFIX="${FFDB_PREFIX_ARG}"


WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       python3 \
       python3-pip \
       python3-setuptools \
       python3-wheel \
       git \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${FFDB_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${FFDB_TAG}" \
  && pip3 install --prefix="${FFDB_PREFIX}" . \
  && add_runtime_dep python3


FROM "${IMAGE}"

ARG FFDB_TAG
ARG FFDB_PREFIX_ARG="/opt/ffdb/${FFDB_TAG}"
ENV FFDB_PREFIX="${FFDB_PREFIX_ARG}"
LABEL ffdb.version="${FFDB_TAG}"

ENV PATH "${FFDB_PREFIX}/bin:${PATH}"
ENV PYTHONPATH "${FFDB_PREFIX}/lib/python3.7/site-packages:${PYTHONPATH}"

COPY --from=ffdb_builder "${FFDB_PREFIX}" "${FFDB_PREFIX}"
COPY --from=ffdb_builder "${APT_REQUIREMENTS_FILE}" /build/apt/ffdb.txt


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
