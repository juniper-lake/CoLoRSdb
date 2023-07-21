FROM python:3.9-buster

LABEL org.opencontainers.image.authors="Juniper Lake <jlake@pacificbiosciences.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		build-essential \
		zlib1g-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
		parallel \
	&& rm -rf /var/lib/apt/lists/*

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts

ARG SOMALIER_VERSION
RUN wget https://github.com/brentp/somalier/releases/download/v${SOMALIER_VERSION}/somalier --directory-prefix /usr/local/bin/ \
  && chmod +x /usr/local/bin/somalier

ARG HTSLIB_VERSION
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
	&& tar -jxvf htslib-${HTSLIB_VERSION}.tar.bz2 --directory /opt \
	&& rm htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt/htslib-${HTSLIB_VERSION} \
	&& make \
	&& make install

ARG LOGURU_VERSION
RUN python3 -m pip install loguru==${LOGURU_VERSION}

ARG PANDAS_VERSION
RUN python3 -m pip install pandas==${PANDAS_VERSION}
