ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y gawk perl \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /opt

ENV SIGNALP_PREFIX="/opt/signalp3"
COPY --from="darcyabjones/signalp3" ${SIGNALP_PREFIX} ${SIGNALP_PREFIX}

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
  && mkdir -p /tmp \
  && sed -i s~/scratch~/tmp~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/bin/perl~$(which perl)~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~AWK=/usr/freeware/bin/gawk~AWK="$(which gawk)"~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/chlorop~${CHLOROP_PREFIX}/chlorop~ ${TARGETP_PREFIX}/targetp \
  && sed -i s~/usr/cbs/bio/bin/signalp~${SIGNALP_PREFIX}/signalp~ ${TARGETP_PREFIX}/targetp


ENV PATH="${TARGETP_PREFIX}:${CHLOROP_PREFIX}:${SIGNALP_PREFIX}:${PATH}"
