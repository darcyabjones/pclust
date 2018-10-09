FROM alpine:3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>

RUN apk add --no-cache perl

WORKDIR /opt
COPY phobius101_linux.tar.gz /opt/phobius101_linux.tar.gz

ENV PHOBIUS_PREFIX="/opt/phobius"

RUN  tar xf phobius101_linux.tar.gz \
  && rm phobius101_linux.tar.gz \
  && mv tmp/tmpTB6e6P/phobius ${PHOBIUS_PREFIX} \
  && ln -sf ${PHOBIUS_PREFIX}/phobius.pl ${PHOBIUS_PREFIX}/phobius

ENV PATH="${PHOBIUS_PREFIX}:${PATH}"
