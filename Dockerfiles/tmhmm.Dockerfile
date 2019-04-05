ARG DEBIAN_VERSION="stretch-20190228-slim"
FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

RUN  apt-get update \
  && apt-get install -y perl \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /opt
COPY tmhmm-2.0c.Linux.tar.gz /opt/tmhmm-2.0c.Linux.tar.gz

ENV TMHMM_PREFIX="/opt/tmhmm"

RUN  tar xf tmhmm-2.0c.Linux.tar.gz \
  && rm tmhmm-2.0c.Linux.tar.gz \
  && mv tmhmm-2.0c ${TMHMM_PREFIX} \
  && sed -i s~/usr/local/bin/perl~$(which perl)~ ${TMHMM_PREFIX}/bin/tmhmm \
  && sed -i s~/usr/local/bin/perl~$(which perl)~ ${TMHMM_PREFIX}/bin/tmhmmformat.pl

ENV PATH="${TMHMM_PREFIX}/bin:${PATH}"
