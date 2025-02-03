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

docker-clean-images:
	docker system prune
