ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

## Config variables
ENV PYTHON3_PREFIX="/opt/python3"
ARG PYTHON3_VERSION="3.7.2"

RUN  apt-get update \
  && apt-get install -y \
       build-essential \
       wget \
       libbz2-dev \
       libffi-dev \
       liblzma-dev \
       libsqlite3-dev \
       libssl1.0-dev \
       zlib1g-dev \
  && rm -rf -- /var/lib/apt/lists/*

WORKDIR /tmp/python3
RUN  wget https://www.python.org/ftp/python/${PYTHON3_VERSION}/Python-${PYTHON3_VERSION}.tgz \
  && tar xfz Python-${PYTHON3_VERSION}.tgz \
  && cd Python-${PYTHON3_VERSION} \
  && ./configure --enable-optimizations --prefix ${PYTHON3_PREFIX} \
  && make \
  && make install


FROM debian:${DEBIAN_VERSION}
LABEL maintainer="darcy.ab.jones@gmail.com"

ENV PYTHON3_PREFIX="/opt/python3"
COPY --from=builder "${PYTHON3_PREFIX}" "${PYTHON3_PREFIX}"

RUN  apt-get update \
  && apt-get install -y \
       libbz2-1.0 \
       libffi6 \
       libsqlite3-0 \
       libssl1.0 \
       lzma \
       zlib1g \
  && rm -rf /var/lib/apt/lists/*

ENV PATH="${PATH}:${PYTHON3_PREFIX}/bin"
ENV LIBRARY_PATH="${LD_LIBRARY_PATH}:${PYTHON3_PREFIX}/lib"
ENV LD_LIBRARY_PATH="${LIBRARY_PATH}"
ENV CPATH="${CPATH}:${PYTHON3_PREFIX}/include"

ENTRYPOINT ["python3"]
CMD ["--help"]
