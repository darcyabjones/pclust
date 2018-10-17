FROM alpine:3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache perl bash

WORKDIR /opt
COPY tmhmm-2.0c.Linux.tar.gz /opt/tmhmm-2.0c.Linux.tar.gz

ENV TMHMM_PREFIX="/opt/tmhmm"

RUN  tar xf tmhmm-2.0c.Linux.tar.gz \
  && rm tmhmm-2.0c.Linux.tar.gz \
  && mv tmhmm-2.0c ${TMHMM_PREFIX} \
  && sed -i s~/usr/local/bin/perl~$(which perl)~ ${TMHMM_PREFIX}/bin/tmhmm \
  && sed -i s~/usr/local/bin/perl~$(which perl)~ ${TMHMM_PREFIX}/bin/tmhmmformat.pl



ENV PATH="${TMHMM_PREFIX}/bin:${PATH}"
