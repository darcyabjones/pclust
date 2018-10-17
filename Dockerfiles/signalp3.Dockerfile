FROM alpine:3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache gawk bash

WORKDIR /opt
COPY signalp-3.0.Linux.tar.Z /opt/signalp-3.0.Linux.tar.Z

ENV SIGNALP_PREFIX="/opt/signalp"

RUN  tar xf signalp-3.0.Linux.tar.Z \
  && rm signalp-3.0.Linux.tar.Z \
  && mv signalp-3.0 ${SIGNALP_PREFIX} \
  && sed -i s~SIGNALP=/usr/opt/signalp-3.0~SIGNALP=${SIGNALP_PREFIX}~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~AWK=/usr/bin/gawk~AWK="$(which gawk)"~ ${SIGNALP_PREFIX}/signalp


ENV PATH="${SIGNALP_PREFIX}:${PATH}"
