.PHONY: *

DOCKER=DOCKER_BUILDKIT=1 docker

all:
	# Nothing to be built

docker-build:
	$(DOCKER) build -t genmod/test --force-rm=true --rm=true -f Dockerfile.pytest .

test: docker-build
	$(DOCKER) run -i -l genmod-test genmod/test -v

test_%: docker-build
	$(DOCKER) run -i -l genmod-test genmod/test -v -k $@

test_dist: docker-build
	$(DOCKER) run -i -l genmod-test --entrypoint=bash genmod/test -c \
		"uv build && \
		pip3 install dist/genmod*.tar.gz && \
		uv run genmod -v annotate --annotate-regions --genome-build 37 examples/test_vcf.vcf"

test_rescore: docker-build
	$(DOCKER) run -i -l genmod-test genmod/test -v -s -o log_cli=true -k test_rescore_with_annotation_suffix 2>&1

test_mivmir: docker-build
	$(DOCKER) run -i -l genmod-test genmod/test -v -s -o log_cli=true -k test_mivmir_minimal_score_config 2>&1

build-export-singularity:
	$(DOCKER) build -t genmod/hasta --force-rm=true --rm=true -f Dockerfile .
	docker save genmod/hasta -o genmod.tar
	singularity build -F genmod.sif docker-archive://genmod.tar
	rm genmod.tar

docker-clean-images:
	docker system prune
