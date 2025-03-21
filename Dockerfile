###############################################
# Dockerfile to build GeMod container image
###############################################
FROM ghcr.io/astral-sh/uv:python3.13-alpine

RUN apk update & apk add build-base zlib-dev
# Copy the project into the image
ADD . /app

# Sync the project into a new virtual environment, using the frozen lockfile
WORKDIR /app
RUN uv sync --frozen

# Add the project venv to PATH to be able to run genmod without changing to the environment
ENV PATH="/app/.venv/bin:$PATH"

ENTRYPOINT ["genmod"]
