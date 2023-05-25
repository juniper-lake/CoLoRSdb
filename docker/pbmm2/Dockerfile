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
	&& rm -rf /var/lib/apt/lists/*

ARG PBMM2_VERSION
RUN wget https://github.com/PacificBiosciences/pbmm2/releases/download/v${PBMM2_VERSION}/pbmm2 --directory-prefix /usr/local/bin/ \
	&& chmod +x /usr/local/bin/pbmm2

ARG PYSAM_VERSION
# Issue with the use_2to3 library required to build pysam for setuptools > 58
RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pysam==${PYSAM_VERSION}

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}":/opt/scripts
