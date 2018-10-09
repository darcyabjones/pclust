FROM alpine:3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache gawk perl

WORKDIR /opt

COPY signalp-3.0.Linux.tar.Z /opt/signalp-3.0.Linux.tar.Z
ENV SIGNALP_PREFIX="/opt/signalp"

RUN  tar xf signalp-3.0.Linux.tar.Z \
  && rm signalp-3.0.Linux.tar.Z \
  && mv signalp-3.0 ${SIGNALP_PREFIX} \
  && sed -i s~SIGNALP=/usr/opt/signalp-3.0~SIGNALP=${SIGNALP_PREFIX}~ ${SIGNALP_PREFIX}/signalp


COPY chlorop-1.1.Linux.tar.Z /opt/chlorop-1.1.Linux.tar.Z
ENV CHLOROP_PREFIX="/opt/chlorop"

RUN  tar xf chlorop-1.1.Linux.tar.Z \
  && rm chlorop-1.1.Linux.tar.Z \
  && mv chlorop-1.1 ${CHLOROP_PREFIX} \
  && sed -i s~/usr/cbs/packages/chlorop/currdist/chlorop-1.1~${CHLOROP_PREFIX}~ ${CHLOROP_PREFIX}/chlorop


COPY targetp-1.1b.Linux.tar.Z /opt/targetp-1.1b.Linux.tar.Z
ENV TARGETP_PREFIX="/opt/targetp"

RUN  tar xf targetp-1.1b.Linux.tar.Z \
  && rm targetp-1.1b.Linux.tar.Z \
  && mv targetp-1.1 ${TARGETP_PREFIX} \
  && sed -i s~/usr/cbs/packages/targetp/currdist/targetp-1.1~${TARGETP_PREFIX}~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/scratch~/tmp~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/bin/perl~perl~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/chlorop~${CHLOROP_PREFIX}/chlorop~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/signalp~${SIGNALP_PREFIX}/signalp~ ${TARGETP_PREFIX}/targetp


ENV PATH="${TARGETP_PREFIX}:${CHLOROP_PREFIX}:${SIGNALP_PREFIX}:${PATH}"
