FROM ubuntu:focal

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
	&& rm -rf /var/lib/apt/lists/*

ARG HIPHASE_VERSION
RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/v${HIPHASE_VERSION}/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
	&& tar -xzvf hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt \
	&& rm hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz \
  && chmod +x /opt/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu/hiphase
	
