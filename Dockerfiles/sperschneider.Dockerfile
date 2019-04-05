ARG DEBIAN_VERSION="stretch-20190228-slim"

FROM debian:${DEBIAN_VERSION} as builder

RUN  apt-get update \
  && apt-get install -y wget unzip \
  && rm -rf -- /var/lib/apt/lists/*

ENV EMBOSS_PREFIX="/opt/EMBOSS"
ENV WEKA_381_PREFIX="/opt/weka-3-8-1"
ENV WEKA_3612_PREFIX="/opt/weka-3-6-12"
ENV EFFECTORP2_PREFIX="/opt/effectorp2"
ENV EFFECTORP_PREFIX="/opt/effectorp"
ENV APOPLASTP_PREFIX="/opt/apoplastp"
ENV LOCALIZER_PREFIX="/opt/localizer"


WORKDIR /opt
RUN  wget http://effectorp.csiro.au/EffectorP_1.0.tar.gz \
  && tar zxf EffectorP_1.0.tar.gz \
  && rm EffectorP_1.0.tar.gz \
  && mv EffectorP_1.0 ${EFFECTORP_PREFIX} \
  && cd ${EFFECTORP_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && unzip -d /opt weka-3-6-12.zip \
  && rm weka-3-6-12.zip \
  && chmod a+x EffectorP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" EffectorP.py \
  && ln -sf "${WEKA_3612_PREFIX}" "${EFFECTORP_PREFIX}/Scripts/weka-3-6-12" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" EffectorP.py \
  && ln -sf ${PWD}/EffectorP.py ${PWD}/EffectorP-1.0.py

WORKDIR /opt
RUN  wget http://localizer.csiro.au/LOCALIZER_1.0.4.tar.gz \
  && tar zxf LOCALIZER_1.0.4.tar.gz \
  && rm LOCALIZER_1.0.4.tar.gz \
  && mv LOCALIZER_1.0.4 ${LOCALIZER_PREFIX} \
  && cd ${LOCALIZER_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && rm weka-3-6-12.zip \
  && chmod a+x LOCALIZER.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" LOCALIZER.py \
  && ln -sf "${WEKA_3612_PREFIX}" "${LOCALIZER_PREFIX}/Scripts/weka-3-6-12" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" LOCALIZER.py

WORKDIR /opt
RUN  wget http://apoplastp.csiro.au/ApoplastP_1.0.1.tar.gz \
  && tar zxf ApoplastP_1.0.1.tar.gz \
  && rm ApoplastP_1.0.1.tar.gz \
  && mv ApoplastP_1.0.1 ${APOPLASTP_PREFIX} \
  && cd ${APOPLASTP_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && unzip -d /opt weka-3-8-1.zip \
  && rm weka-3-8-1.zip \
  && chmod a+x ApoplastP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" ApoplastP.py \
  && ln -sf "${WEKA_381_PREFIX}" "${APOPLASTP_PREFIX}/Scripts/weka-3-8-1" \
  && sed -i "s~SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'~'${EMBOSS_PREFIX}/bin/'~" ApoplastP.py

WORKDIR /opt
RUN  wget http://effectorp.csiro.au/EffectorP_2.0.tar.gz \
  && tar zxf EffectorP_2.0.tar.gz \
  && rm EffectorP_2.0.tar.gz \
  && mv EffectorP_2.0 ${EFFECTORP2_PREFIX} \
  && cd ${EFFECTORP2_PREFIX}/Scripts \
  && rm emboss-latest.tar.gz \
  && rm weka-3-8-1.zip \
  && chmod a+x EffectorP.py \
  && sed -i "s~/usr/bin/python~/usr/bin/env python~" EffectorP.py \
  && ln -sf "${WEKA_381_PREFIX}" "${EFFECTORP2_PREFIX}/Scripts/weka-3-8-1" \
  && mkdir "${EFFECTORP2_PREFIX}/Scripts/EMBOSS-6.5.7" \
  && ln -sf "${EMBOSS_PREFIX}/bin" "${EFFECTORP2_PREFIX}/Scripts/EMBOSS-6.5.7/emboss" \
  && ln -sf "${PWD}/EffectorP.py" "${PWD}/EffectorP-2.0.py"


FROM debian:${DEBIAN_VERSION}

LABEL maintainer="darcy.ab.jones@gmail.com"

ENV EMBOSS_PREFIX="/opt/EMBOSS"
ENV WEKA_381_PREFIX="/opt/weka-3-8-1"
ENV WEKA_3612_PREFIX="/opt/weka-3-6-12"
ENV EFFECTORP2_PREFIX="/opt/effectorp2"
ENV EFFECTORP_PREFIX="/opt/effectorp"
ENV APOPLASTP_PREFIX="/opt/apoplastp"
ENV LOCALIZER_PREFIX="/opt/localizer"

COPY --from="darcyabjones/emboss" ${EMBOSS_PREFIX} ${EMBOSS_PREFIX}
COPY --from=builder ${WEKA_3612_PREFIX} ${WEKA_3612_PREFIX}
COPY --from=builder ${WEKA_381_PREFIX} ${WEKA_381_PREFIX}
COPY --from=builder ${EFFECTORP2_PREFIX} ${EFFECTORP2_PREFIX}
COPY --from=builder ${EFFECTORP_PREFIX} ${EFFECTORP_PREFIX}
COPY --from=builder ${APOPLASTP_PREFIX} ${APOPLASTP_PREFIX}
COPY --from=builder ${LOCALIZER_PREFIX} ${LOCALIZER_PREFIX}

# The man thing is required to get java installed on the minimal distro
RUN  if [ ! -d /usr/share/man/man1 ]; then \
       mkdir -p /usr/share/man/man1; \
     fi; \
     apt-get update \
     && apt-get install -y \
       default-jre \
       python2.7 \
       python-numpy \
       python-biopython \
  && rm -rf -- /var/lib/apt/lists/*

ENV PATH="${EMBOSS_PREFIX}/bin:${EFFECTORP2_PREFIX}/Scripts:${EFFECTORP_PREFIX}/Scripts:${APOPLASTP_PREFIX}/Scripts:${LOCALIZER_PREFIX}/Scripts:${PATH}"
