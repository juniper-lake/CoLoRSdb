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
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
		liblzma-dev \
		libbz2-dev \
		libcurl4-openssl-dev \
		parallel \
	&& rm -rf /var/lib/apt/lists/*

ARG PBSV_VERSION
RUN wget https://github.com/PacificBiosciences/pbsv/releases/download/v${PBSV_VERSION}/pbsv --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pbsv

ARG SAMTOOLS_VERSION
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 --directory /opt \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt/samtools-${SAMTOOLS_VERSION} \
	&& make \
	&& make install
