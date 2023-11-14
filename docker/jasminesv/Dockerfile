FROM continuumio/miniconda3:4.12.0

LABEL org.opencontainers.image.authors="Juniper Lake <jlake@pacificbiosciences.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

# Create environment and add to path
ARG JASMINESV_VERSION
ARG JASMINESV_BUILD
RUN conda create --yes --channel bioconda --name jasminesv jasminesv=${JASMINESV_VERSION}=${JASMINESV_BUILD}
ENV PATH /opt/conda/envs/${IMAGE_NAME}/bin:$PATH

# Add to ~/.bashrc in case it needs to be sourced
RUN echo "source activate ${IMAGE_NAME}" >> ~/.bashrc
