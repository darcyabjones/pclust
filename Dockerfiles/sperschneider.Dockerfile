FROM python:2.7.15-alpine3.8
MAINTAINER Darcy Jones <darcy.ab.jones@gmail.com>
RUN apk add --no-cache build-base bash openjdk7="7.181.2.6.14-r0" perl \
  && pip install --no-cache-dir numpy biopython

ENV EMBOSS_PREFIX="/opt/EMBOSS-6.5.7"
ENV WEKA_381_PREFIX="/opt/weka-3-8-1"
ENV WEKA_3612_PREFIX="/opt/weka-3-6-12"
ENV EFFECTORP_PREFIX="/opt/effectorp"
ENV APOPLASTP_PREFIX="/opt/apoplastp"
ENV LOCALIZER_PREFIX="/opt/localizer"

WORKDIR /opt
RUN  wget http://effectorp.csiro.au/EffectorP_2.0.tar.gz \
  && tar zxf EffectorP_2.0.tar.gz \
  && rm EffectorP_2.0.tar.gz \
  && mv EffectorP_2.0 ${EFFECTORP_PREFIX} \
  && cd ${EFFECTORP_PREFIX}/Scripts \
  && mv ${EFFECTORP_PREFIX}/Scripts/emboss-latest.tar.gz /opt \
  && mv ${EFFECTORP_PREFIX}/Scripts/weka-3-8-1.zip /opt \
  && chmod a+x EffectorP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" EffectorP.py \
  \
  && cd /opt \
  && wget http://apoplastp.csiro.au/ApoplastP_1.0.1.tar.gz \
  && tar zxf ApoplastP_1.0.1.tar.gz \
  && rm ApoplastP_1.0.1.tar.gz \
  && mv ApoplastP_1.0.1 ${APOPLASTP_PREFIX} \
  && cd ${APOPLASTP_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && rm weka-3-8-1.zip \
  && chmod a+x ApoplastP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" ApoplastP.py \
  \
  && cd /opt \
  && wget http://localizer.csiro.au/LOCALIZER_1.0.4.tar.gz \
  && tar zxf LOCALIZER_1.0.4.tar.gz \
  && rm LOCALIZER_1.0.4.tar.gz \
  && mv LOCALIZER_1.0.4 ${LOCALIZER_PREFIX} \
  && cd ${LOCALIZER_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && mv weka-3-6-12.zip /opt \
  && chmod a+x LOCALIZER.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" LOCALIZER.py \
  \
  && cd /opt \
  && tar zxf emboss-latest.tar.gz \
  && cd EMBOSS-6.5.7 \
  && ./configure --without-x \
  && make \
  && cd /opt \
  && unzip weka-3-8-1.zip \
  && unzip weka-3-6-12.zip \
  && rm weka-3-6-12.zip \
  && rm weka-3-8-1.zip \
  && rm emboss-latest.tar.gz \
  && ln -sf ${EMBOSS_PREFIX} ${EFFECTORP_PREFIX}/Scripts \
  && ln -sf ${EMBOSS_PREFIX} ${APOPLASTP_PREFIX}/Scripts \
  && ln -sf ${EMBOSS_PREFIX} ${LOCALIZER_PREFIX}/Scripts \
  && ln -sf ${WEKA_381_PREFIX} ${EFFECTORP_PREFIX}/Scripts \
  && ln -sf ${WEKA_381_PREFIX} ${APOPLASTP_PREFIX}/Scripts \
  && ln -sf ${WEKA_3612_PREFIX} ${LOCALIZER_PREFIX}/Scripts


ENV PATH="${EFFECTORP_PREFIX}/Scripts:${APOPLASTP_PREFIX}/Scripts:${LOCALIZER_PREFIX}/Scripts:${PATH}"
