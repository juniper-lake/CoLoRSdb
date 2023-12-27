# Dockers

This workflow requires docker images defined in [the docker directory](docker) that are hosted in the [CoLoRSdb quay.io repo](https://quay.io/colorsdb).

## Directory structure

Docker image definitions found in [the docker directory](docker) define a [`build.env` file](#the-buildenv-file) and a Dockerfile that gets populated with variables from the [`build.env` file](#the-buildenv-file).

Example directory structure:
```
docker/
├── bcftools/
│   ├── build.env
│   └── Dockerfile
└──somalier/
    ├── build.env
    ├── Dockerfile
    └── scripts/

```

## The `build.env` file

Each target image is defined via the presence of a `build.env` file, which is used to specify the name and version tag for the corresponding Docker image. It must contain at minimum the following variables:

- `IMAGE_NAME`: specifies the name of the built image
- `IMAGE_TAG`: specifies the tag of the built image

All variables defined in the `build.env` file will be made available as build arguments during Docker image build.

The `IMAGE_TAG` variable can be built using other variables defined in the `build.env` file, as long as those other variables are defined before `IMAGE_TAG`. For example, the following `IMAGE_TAG` would be set to `0.7.8_1.15`:

```
# Tool versions
BWA_VERSION=0.7.8
SAMTOOLS_VERSION=1.15

# Image info
IMAGE_NAME=bwa_samtools
IMAGE_TAG=${BWA_VERSION}_${SAMTOOLS_VERSION}
```

The `Dockerfile` corresponding to a `build.env` must either:
1. Be named `Dockerfile` and exist in the same directory as the `build.env` file, or
2. Be specified using the `DOCKERFILE` variable in the `build.env` file (e.g. for image builds that reuse a common Dockerfile). The path to the Dockerfile should be relative to the directory containing the `build.env` file.

## Building Docker images

Docker images can be built using the [build_docker_images](util/build_docker_images) utility script.


### Build a single image

```bash
./util/build_docker_images -d docker/smrtcell_stats
```

#### Build all images in the `docker` directory

```bash
./util/build_docker_images -d docker
```

#### Build and push all images in the `docker` directory

```bash
./util/build_docker_images -d docker -p
```

#### Build and push all images in the `docker` directory, using the `colorsdb` container registry

```bash
./util/build_docker_images -d docker -p -c colorsdb
```

## Tool versions and Docker images

The following docker images are hosted in the [CoLoRSdb quay.io repo](https://quay.io/colorsdb).

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Docker images used in the workflow are pegged to specific versions by referring to their digests rather than tags. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`.

| Image | Major tool versions | Links |
| :- | :- | :- |
| bcftools | <ul><li>[bcftools 1.14](https://github.com/samtools/bcftools/releases/tag/1.14)</li><li>[htslib 1.18](https://github.com/samtools/htslib/releases/tag/1.18)</li></ul> | [Dockerfile](docker/bcftools) |
| deepvariant | <ul><li>[deepvariant 1.5.0](https://github.com/google/deepvariant/releases/tag/v1.5.0)</li></ul> | [DeepVariant GitHub](https://github.com/google/deepvariant) |
| glnexus | <ul><li>[glnexus v1.4.1](https://github.com/dnanexus-rnd/GLnexus/releases/tag/v1.4.1)</li></ul> | [GLnexus GitHub](https://github.com/dnanexus-rnd/GLnexus) |
| hificnv | <ul><li>[HiFiCNV v0.1.6](https://github.com/PacificBiosciences/HiFiCNV/releases/tag/v0.1.6)</li><li>[bcftools 1.17](https://github.com/samtools/bcftools/releases/tag/1.17)</li><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li></ul> | [Dockerfile](docker/hificnv) |
| htslib | <ul><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li></ul> | [Dockerfile](docker/htslib) |
| jasminesv | <ul><li>[jasmine 1.1.5](https://github.com/mkirsche/Jasmine/releases/tag/v1.1.5)</li></ul> | [Dockerfile](docker/jasminesv) |
| mosdepth | <ul><li>[mosdepth 0.2.9](https://github.com/brentp/mosdepth/releases/tag/v0.2.9)</li></ul> | [Dockerfile](docker/mosdepth) |
| pbmm2 | <ul><li>[pbmm2 1.10.0](https://github.com/PacificBiosciences/pbmm2/releases/tag/v1.10.0)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li></ul> | [Dockerfile](docker/pbmm2) |
| pbsv | <ul><li>[pbsv 2.9.0](https://github.com/PacificBiosciences/pbsv/releases/tag/v2.9.0)</li></ul> | [Dockerfile](docker/pbsv) |
| peddy | <ul><li>[peddy 0.4.8](https://github.com/brentp/peddy/releases/tag/v0.4.8)</li> | [Dockerfile](docker/peddy)
| samtools | <ul><li>[samtools 1.17](https://github.com/samtools/samtools/releases/tag/1.17)</li></ul> | [Dockerfile](docker/samtools) |
| sniffles | <ul><li>[sniffles 2.2](https://github.com/fritzsedlazeck/Sniffles/releases/tag/v2.2)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li></ul> | [Dockerfile](docker/sniffles)
| somalier | <ul><li>[somalier 0.2.16](https://github.com/brentp/somalier/releases/tag/v0.2.16)</li><li>python 3.9; custom scripts</li><li>[loguru 0.7.0](https://github.com/Delgan/loguru/releases/tag/0.7.0)</li></ul> | [Dockerfile](docker/somalier)
| trgt | <ul><li>[trgt 0.5.0](https://github.com/PacificBiosciences/trgt/releases/tag/v0.5.0)</li><li>[bcftools 1.16](https://github.com/samtools/bcftools/releases/tag/1.16)</li></ul> | [Dockerfile](docker/trgt) |
| vcfparser | <ul><li>python 3.9; custom scripts</li><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li><li>[loguru 0.7.0](https://github.com/Delgan/loguru/releases/tag/0.7.0)</li><li>[bedtools 2.27.1](https://github.com/arq5x/bedtools2/releases/tag/v2.27.1)</li></ul> | [Dockerfile](docker/vcfparser) |
