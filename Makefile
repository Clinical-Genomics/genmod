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

docker-clean-images:
	docker system prune
