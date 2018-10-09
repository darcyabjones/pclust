FROM alpine:3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache perl

WORKDIR /opt
COPY signalp-4.1f.Linux.tar.gz /opt/signalp-4.1f.Linux.tar.gz

ENV SIGNALP_PREFIX="/opt/signalp"

RUN  tar xf signalp-4.1f.Linux.tar.gz \
  && rm signalp-4.1f.Linux.tar.gz \
  && mv signalp-4.1 ${SIGNALP_PREFIX} \
  && sed -i s~/usr/cbs/bio/src/signalp-4.1~${SIGNALP_PREFIX}~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~/var/tmp~/tmp~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~MAX_ALLOWED_ENTRIES=10000~MAX_ALLOWED_ENTRIES=999999999999999~ ${SIGNALP_PREFIX}/signalp


ENV PATH="${SIGNALP_PREFIX}:${PATH}"
