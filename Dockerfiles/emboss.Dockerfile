ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

ENV EMBOSS_PREFIX="/opt/EMBOSS"
ARG EMBOSS_VERSION="6.5.7"
ARG EMBOSS_URL="ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz"

RUN  apt-get update \
  && apt-get install -y build-essential wget \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/emboss
RUN  wget ${EMBOSS_URL} \
  && tar xfz EMBOSS-${EMBOSS_VERSION}.tar.gz \
  && cd EMBOSS-${EMBOSS_VERSION} \
  && ./configure --without-x --prefix ${EMBOSS_PREFIX} \
  && make \
  && make install


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV EMBOSS_PREFIX="/opt/EMBOSS"
COPY --from=builder ${EMBOSS_PREFIX} ${EMBOSS_PREFIX}

ENV PATH="${PATH}:${EMBOSS_PREFIX}/bin"
