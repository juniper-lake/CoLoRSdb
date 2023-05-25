FROM ubuntu:focal

LABEL org.opencontainers.image.authors="Juniper Lake <jlake@pacificbiosciences.com>"
LABEL org.opencontainers.image.url="https://github.com/juniper-lake/CoLoRS"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		parallel \
	&& rm -rf /var/lib/apt/lists/*

ARG DATAMASH_VERSION
RUN wget https://ftp.gnu.org/gnu/datamash/datamash-${DATAMASH_VERSION}.tar.gz \
	&& tar -zxvf datamash-${DATAMASH_VERSION}.tar.gz --directory /opt \
	&& rm datamash-${DATAMASH_VERSION}.tar.gz
RUN cd /opt/datamash-${DATAMASH_VERSION} \
	&& ./configure \
	&& make \
	&& make check \
	&& make install
  