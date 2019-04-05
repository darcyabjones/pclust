ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y perl \
  && rm -rf -- /var/lib/apt/lists/*

ENV SIGNALP_PREFIX="/opt/signalp4"

WORKDIR /opt
COPY signalp-4.1f.Linux.tar.gz /opt/signalp-4.1f.Linux.tar.gz

RUN  tar xf signalp-4.1f.Linux.tar.gz \
  && rm signalp-4.1f.Linux.tar.gz \
  && mv signalp-4.1 ${SIGNALP_PREFIX} \
  && sed -i s~/usr/cbs/bio/src/signalp-4.1~${SIGNALP_PREFIX}~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~/var/tmp~/tmp~ ${SIGNALP_PREFIX}/signalp \
  && sed -i s~MAX_ALLOWED_ENTRIES=10000~MAX_ALLOWED_ENTRIES=999999999999999~ ${SIGNALP_PREFIX}/signalp \
  && ln -sf ${SIGNALP_PREFIX}/signalp ${SIGNALP_PREFIX}/signalp-4.1


ENV PATH="${SIGNALP_PREFIX}:${PATH}"
