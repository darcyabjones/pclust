ARG IMAGE
ARG MMSEQS_IMAGE

FROM "${MMSEQS_IMAGE}" as mmseqs_builder

FROM "${IMAGE}" as fasttree_builder

ARG FASTTREE_VERSION="2.1.11"
ARG FASTTREE_URL="http://www.microbesonline.org/fasttree/FastTree-2.1.11.c"
ARG FASTTREE_PREFIX_ARG="/opt/fasttree/${FASTTREE_VERSION}"
ENV FASTTREE_PREFIX="${FASTTREE_PREFIX_ARG}"

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
  && wget -O FastTree.c "${FASTTREE_URL}" \
  && gcc -msse4 -mno-avx -O3 -finline-functions -funroll-loops -Wall -o FastTree-sse4 FastTree.c -lm \
  && gcc -msse4 -mno-avx -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP-sse4 FastTree.c -lm \
  && gcc -mavx2 -O3 -finline-functions -funroll-loops -Wall -o FastTree-avx2 FastTree.c -lm \
  && gcc -mavx2 -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP-avx2 FastTree.c -lm \
  && mkdir -p "${FASTTREE_PREFIX}/bin" \
  && mv FastTree-* "${FASTTREE_PREFIX}/bin" \
  && mv FastTreeMP-* "${FASTTREE_PREFIX}/bin" \
  && echo '#!/usr/bin/env bash'                                  > "${FASTTREE_PREFIX}/bin/FastTree" \
  && echo 'if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then' >> "${FASTTREE_PREFIX}/bin/FastTree" \
  && echo '    exec "${FASTTREE_PREFIX}/bin/FastTree-avx2" "$@"'    >> "${FASTTREE_PREFIX}/bin/FastTree" \
  && echo 'else'                                                >> "${FASTTREE_PREFIX}/bin/FastTree" \
  && echo '    exec "${FASTTREE_PREFIX}/bin/FastTree-sse4" "$@"'    >> "${FASTTREE_PREFIX}/bin/FastTree" \
  && echo 'fi'                                                  >> "${FASTTREE_PREFIX}/bin/FastTree" \
  && sed 's/FastTree-/FastTreeMP-/g' "${FASTTREE_PREFIX}/bin/FastTree" >"${FASTTREE_PREFIX}/bin/FastTreeMP" \
  && chmod a+x "${FASTTREE_PREFIX}/bin/FastTree" \
  && chmod a+x "${FASTTREE_PREFIX}/bin/FastTreeMP" \
  && add_runtime_dep \
       gawk \
       bash \
       grep \
       libgomp1


FROM "${IMAGE}"

ARG FASTTREE_VERSION="2.1.11"
ARG FASTTREE_PREFIX_ARG="/opt/fasttree/${FASTTREE_VERSION}"
ENV FASTTREE_PREFIX="${FASTTREE_PREFIX_ARG}"
LABEL fasttree.version="${FASTTREE_VERSION}"

COPY --from=fasttree_builder "${FASTTREE_PREFIX}" "${FASTTREE_PREFIX}"
COPY --from=fasttree_builder "${APT_REQUIREMENTS_FILE}" /build/apt/fasttree.txt

ENV PATH="${FASTTREE_PREFIX}/bin:${PATH}"


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
