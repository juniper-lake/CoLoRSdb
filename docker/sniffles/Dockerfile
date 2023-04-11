FROM python:3.9-buster

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

ARG SNIFFLES_VERSION
RUN python3 -m pip install paraphase==${SNIFFLES_VERSION}

ARG PYSAM_VERSION
RUN python3 -m pip install setuptools==58
RUN python3 -m pip install pysam==${PYSAM_VERSION}
