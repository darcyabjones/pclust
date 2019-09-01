PWD=$(shell pwd)
DEBIAN_VERSION=buster-20190708-slim
TAG=pclust-v0.0.1
MAINTAINER=darcy.ab.jones@gmail.com
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs

MMSEQS_TAG=9-d36de
MMSEQS_REPO=https://github.com/soedinglab/MMseqs2.git
MMSEQS_PREFIX_ARG=$(INSTALL_BASE)/mmseqs/$(MMSEQS_TAG)

HHSUITE_TAG=v3.2.0
HHSUITE_REPO=https://github.com/soedinglab/hh-suite.git
HHSUITE_PREFIX_ARG=$(INSTALL_BASE)/hhsuite/$(HHSUITE_TAG)
HHSUITE_CMAKE_OPTIONS=

MAFFT_VERSION=7.407
MAFFT_URL=https://mafft.cbrc.jp/alignment/software/mafft-$(MAFFT_VERSION)-without-extensions-src.tgz
MAFFT_PREFIX_ARG=$(INSTALL_BASE)/mafft/$(MAFFT_VERSION)

FFDB_TAG=v0.0.1
FFDB_REPO=https://github.com/darcyabjones/ffdb.git
FFDB_PREFIX_ARG=$(INSTALL_BASE)/ffdb/$(FFDB_TAG)


GENERAL=mmseqs python3 posix hhsuite mafft
OPEN=$(GENERAL)

all: docker/all singularity/all

# Docker targets
docker: docker/panann
docker/all: docker/general

docker/general: $(addprefix docker/,$(GENERAL))

# Singularity targets
singularity: singularity/panann.sif
singularity/all: singularity/general

singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))

# Pull containers for dockerhub
dockerpull: dockerpull/panann
dockerpull/all: dockerpull/general
dockerpull/general: $(addprefix dockerpull/,$(GENERAL))


# Utils
clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean -a
	rm -rf -- tmp
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	singularity cache clean -a
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --build-arg MAINTAINER="$(MAINTAINER)" \
	  --build-arg PIPELINE_VERSION="$(TAG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/posix: posix.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mmseqs: mmseqs.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_REPO="$(MMSEQS_REPO)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/hhsuite: hhsuite.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HHSUITE_TAG="$(HHSUITE_TAG)" \
	  --build-arg HHSUITE_REPO="$(HHSUITE_REPO)" \
	  --build-arg HHSUITE_PREFIX_ARG="$(HHSUITE_PREFIX_ARG)" \
	  --build-arg HHSUITE_CMAKE_OPTIONS="$(HHSUITE_CMAKE_OPTIONS)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/ffdb: ffdb.Dockerfile docker/base docker/hhsuite
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FFDB_TAG="$(FFDB_TAG)" \
	  --build-arg FFDB_REPO="$(FFDB_REPO)" \
	  --build-arg FFDB_PREFIX_ARG="$(FFDB_PREFIX_ARG)" \
	  --build-arg HHSUITE_IMAGE="$(call make_image_name,hhsuite)" \
	  --build-arg HHSUITE_TAG="$(HHSUITE_TAG)" \
	  --build-arg HHSUITE_PREFIX_ARG="$(HHSUITE_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mafft: mafft.Dockerfile docker/base docker/mmseqs
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MAFFT_VERSION="$(MAFFT_VERSION)" \
	  --build-arg MAFFT_PREFIX_ARG="$(MAFFT_PREFIX_ARG)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@
