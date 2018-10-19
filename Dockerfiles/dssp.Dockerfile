FROM alpine:3.8 as builder
RUN apk add --no-cache build-base boost-dev automake autoconf bzip2-dev zlib-dev

WORKDIR /opt
RUN    wget https://github.com/cmbi/xssp/archive/3.0.5.tar.gz \
    && tar -zxf 3.0.5.tar.gz \
    && rm 3.0.5.tar.gz \
    && cd xssp-3.0.5 \
    && ./autogen.sh \
    && ./configure --prefix /opt/xssp \
    && sed -i '30a #ifndef uint\n#define uint unsigned int\n#endif' src/progress.cpp \
    && make \
    && make install \
    && rm -rf -- xssp-3.0.5

FROM alpine:3.8

RUN apk add --no-cache bash libstdc++ libbz2 boost-iostreams boost-program_options boost-thread

COPY --from=builder /opt/xssp /opt/xssp
ENV PATH=/opt/xssp/bin:${PATH}  
