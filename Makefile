.PHONY: *

DOCKER=DOCKER_BUILDKIT=1 docker

all:
	# Nothing to be built

docker-build:
	$(DOCKER) build -t genmod/test --force-rm=true --rm=true -f Dockerfile.pytest .

test: docker-build
	$(DOCKER) run -it -l genmod-test genmod/test

test_%: docker-build
	$(DOCKER) run -it -l genmod-test genmod/test -k $@

test_dist: docker-build
	$(DOCKER) run -i -l genmod-test --entrypoint=bash genmod/test -c \
		"python3 -m pip install --upgrade pip && \
		pip3 install build && \
		python3 -m build && \
		python3 -m venv genmodinstall && \
		. genmodinstall/bin/activate && \
		pip3 install dist/genmod*.tar.gz && \
		genmod -v annotate --annotate-regions --genome-build 37 examples/test_vcf.vcf"

docker-clean-images:
	docker system prune
