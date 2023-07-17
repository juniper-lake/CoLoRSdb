FROM ubuntu:focal

LABEL org.opencontainers.image.authors="Juniper Lake <jlake@pacificbiosciences.com>"

RUN apt-get -qq update \
    && apt-get -qq install \
        wget \
        build-essential \
        zlib1g-dev \
        libncurses5-dev \
        liblzma-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        python3 \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

ARG PEDDY_VERSION
RUN wget https://github.com/brentp/peddy/archive/refs/tags/v${PEDDY_VERSION}.tar.gz \
    && tar -xvf v${PEDDY_VERSION}.tar.gz --directory /opt \
    && rm v${PEDDY_VERSION}.tar.gz

RUN cd /opt/peddy-${PEDDY_VERSION} \
    && python3 -m pip install -r requirements.txt \
    && python3 -m pip install --editable .
