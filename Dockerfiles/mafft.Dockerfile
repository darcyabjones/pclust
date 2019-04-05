ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION} as builder

## Config variables
ENV MAFFT_PREFIX="/opt/mafft"
ARG MAFFT_VERSION="7.407"

RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/mafft
RUN  wget https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
  && tar zxf mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
  && cd mafft-${MAFFT_VERSION}-without-extensions \
  && sed -i "s~PREFIX = /usr/local~PREFIX = ${MAFFT_PREFIX}~" core/Makefile \
  && cd core \
  && make clean \
  && make install


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV MAFFT_PREFIX="/opt/mafft"

COPY --from=builder ${MAFFT_PREFIX} ${MAFFT_PREFIX}
ENV PATH="${MAFFT_PREFIX}/bin:${PATH}"

ENTRYPOINT ["mafft"]
CMD ["--help"]
