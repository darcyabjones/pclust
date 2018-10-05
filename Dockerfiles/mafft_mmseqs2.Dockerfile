from alpine:3.8 as mafftbuilder
RUN apk add --no-cache build-base

ENV MAFFT_PREFIX="/opt/mafft"
ARG MAFFT_VERSION="7.407"

# Download, build, and install MPICH
WORKDIR /tmp/mafft
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
      && tar zxf mafft-${MAFFT_VERSION}-without-extensions-src.tgz \
      && cd mafft-${MAFFT_VERSION}-without-extensions \
      && sed -i "s~PREFIX = /usr/local~PREFIX = ${MAFFT_PREFIX}~" core/Makefile \
      && cd core \
      && make clean \
      && make install


FROM soedinglab/mmseqs2:version-5

MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache bash libstdc++ libgomp

# Include HHblits
COPY --from=mafftbuilder /opt/mafft /opt/mafft
ENV PATH="/opt/mafft/bin:${PATH}"
