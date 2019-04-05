ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y gawk \
  && rm -rf -- /var/lib/apt/lists/*

ENV SIGNALP_PREFIX="/opt/signalp3"

WORKDIR /opt
COPY signalp-3.0.Linux.tar.Z /opt/signalp-3.0.Linux.tar.Z

RUN  tar xf signalp-3.0.Linux.tar.Z \
  && rm signalp-3.0.Linux.tar.Z \
  && mv signalp-3.0 ${SIGNALP_PREFIX} \
  && sed -i s~SIGNALP=/usr/opt/signalp-3.0~SIGNALP=${SIGNALP_PREFIX}~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~AWK=/usr/bin/gawk~AWK="$(which gawk)"~ ${SIGNALP_PREFIX}/signalp \
  && ln -sf ${SIGNALP_PREFIX}/signalp ${SIGNALP_PREFIX}/signalp-3


ENV PATH="${SIGNALP_PREFIX}:${PATH}"
