PWD=$(shell pwd)
DEBIAN_VERSION=buster-20190812-slim
PIPELINE_NAME=pclust
PIPELINE_VERSION=v0.0.3
TAG=$(PIPELINE_NAME)-$(PIPELINE_VERSION)
MAINTAINER=darcy.ab.jones@gmail.com
SINGULARITY_USER=darcyabjones
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs

MMSEQS_TAG=10-6d92c
MMSEQS_REPO=https://github.com/soedinglab/MMseqs2.git
MMSEQS_PREFIX_ARG=$(INSTALL_BASE)/mmseqs/$(MMSEQS_TAG)

HHSUITE_TAG=v3.2.0
HHSUITE_REPO=https://github.com/soedinglab/hh-suite.git
HHSUITE_PREFIX_ARG=$(INSTALL_BASE)/hhsuite/$(HHSUITE_TAG)
HHSUITE_CMAKE_OPTIONS=

MAFFT_VERSION=7.407
MAFFT_URL=https://mafft.cbrc.jp/alignment/software/mafft-$(MAFFT_VERSION)-without-extensions-src.tgz
MAFFT_PREFIX_ARG=$(INSTALL_BASE)/mafft/$(MAFFT_VERSION)

FASTTREE_VERSION=2.1.11
FASTTREE_URL=http://www.microbesonline.org/fasttree/FastTree-2.1.11.c
FASTTREE_PREFIX_ARG=$(INSTALL_BASE)/fasttree/$(FASTTREE_VERSION)

FFDB_TAG=93a4a72

FFDB_REPO=https://github.com/darcyabjones/ffdb.git
FFDB_PREFIX_ARG=$(INSTALL_BASE)/ffdb/$(FFDB_TAG)


GENERAL=base mmseqs python3 hhsuite mafft fasttree ffdb
OPEN=$(GENERAL)

all: docker/all singularity/all

# Docker targets
docker: docker/pclust
	touch $@

docker/all: docker/general
	touch $@

docker/general: $(addprefix docker/,$(GENERAL))
	touch $@

# Singularity targets
singularity: singularity/pclust.sif
	touch $@

singularity/all: singularity/general
	touch $@

singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))
	touch $@

# Pull containers for dockerhub
dockerpull: dockerpull/pclust
	touch docker/$(notdir $@)

dockerpull/all: dockerpull/general
	touch docker/$(notdir $@)

dockerpull/general: $(addprefix dockerpull/,$(GENERAL))
	touch docker/$(notdir $@)

# Utils
clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean --force
	rm -rf -- tmp
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	singularity cache clean --force
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	@mkdir -p "${PWD}/tmp"
	@mkdir -p "${PWD}/tmp/cache"
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

singularity/push: singularity/pclust.sif
	singularity sign $<
	singularity push $< library://$(SINGULARITY_USER)/default/$(PIPELINE_NAME):$(PIPELINE_VERSION)

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/push: docker/pclust
	docker login --username $(DOCKER_USER)


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


docker/ffdb: ffdb.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FFDB_TAG="$(FFDB_TAG)" \
	  --build-arg FFDB_REPO="$(FFDB_REPO)" \
	  --build-arg FFDB_PREFIX_ARG="$(FFDB_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mmseqs: mmseqs.Dockerfile docker/base docker/ffdb
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_REPO="$(MMSEQS_REPO)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg FFDB_IMAGE="$(call make_image_name,ffdb)" \
	  --build-arg FFDB_TAG="$(FFDB_TAG)" \
	  --build-arg FFDB_PREFIX_ARG="$(FFDB_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/hhsuite: hhsuite.Dockerfile docker/base docker/ffdb
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HHSUITE_TAG="$(HHSUITE_TAG)" \
	  --build-arg HHSUITE_REPO="$(HHSUITE_REPO)" \
	  --build-arg HHSUITE_PREFIX_ARG="$(HHSUITE_PREFIX_ARG)" \
	  --build-arg HHSUITE_CMAKE_OPTIONS="$(HHSUITE_CMAKE_OPTIONS)" \
	  --build-arg FFDB_IMAGE="$(call make_image_name,ffdb)" \
	  --build-arg FFDB_TAG="$(FFDB_TAG)" \
	  --build-arg FFDB_PREFIX_ARG="$(FFDB_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/mafft: mafft.Dockerfile docker/base docker/mmseqs
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MAFFT_VERSION="$(MAFFT_VERSION)" \
	  --build-arg MAFFT_URL="$(MAFFT_URL)" \
	  --build-arg MAFFT_PREFIX_ARG="$(MAFFT_PREFIX_ARG)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/fasttree: fasttree.Dockerfile docker/base docker/mmseqs
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg FASTTREE_VERSION="$(FASTTREE_VERSION)" \
	  --build-arg FASTTREE_URL="$(FASTTREE_URL)" \
	  --build-arg FASTTREE_PREFIX_ARG="$(FASTTREE_PREFIX_ARG)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/pclust: pclust.Dockerfile docker/all
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg PYTHON3_IMAGE="$(call make_image_name,python3)" \
	  --build-arg MMSEQS_IMAGE="$(call make_image_name,mmseqs)" \
	  --build-arg MMSEQS_TAG="$(MMSEQS_TAG)" \
	  --build-arg MMSEQS_PREFIX_ARG="$(MMSEQS_PREFIX_ARG)" \
	  --build-arg HHSUITE_IMAGE="$(call make_image_name,hhsuite)" \
	  --build-arg HHSUITE_TAG="$(HHSUITE_TAG)" \
	  --build-arg HHSUITE_PREFIX_ARG="$(HHSUITE_PREFIX_ARG)" \
	  --build-arg MAFFT_IMAGE="$(call make_image_name,mafft)" \
	  --build-arg MAFFT_VERSION="$(MAFFT_VERSION)" \
	  --build-arg MAFFT_PREFIX_ARG="$(MAFFT_PREFIX_ARG)" \
	  --build-arg FASTTREE_IMAGE="$(call make_image_name,fasttree)" \
	  --build-arg FASTTREE_VERSION="$(FASTTREE_VERSION)" \
	  --build-arg FASTTREE_PREFIX_ARG="$(FASTTREE_PREFIX_ARG)" \
	  --build-arg FFDB_IMAGE="$(call make_image_name,ffdb)" \
	  --build-arg FFDB_TAG="$(FFDB_TAG)" \
	  --build-arg FFDB_PREFIX_ARG="$(FFDB_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@
