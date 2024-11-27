###############################################
# Dockerfile to build GeMod container image
###############################################
#FROM ghcr.io/astral-sh/uv:python3.13-alpine
FROM ghcr.io/astral-sh/uv:python3.13-bookworm

# Copy the project into the image
ADD . /app

# Sync the project into a new environment, using the frozen lockfile
WORKDIR /app
RUN uv sync --frozen

ENTRYPOINT ["uv", "run", "genmod"]
