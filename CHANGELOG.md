# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

About changelog [here](https://keepachangelog.com/en/1.0.0/)

Please add a new candidate release at the top after changing the latest one. Feel free to copy paste from the "squash and commit" box that gets generated when creating PRs

Try to use the following format:

## [unreleased]
### Fixed
- The optional fields Source and Version are allowed in the VCF header([#106](https://github.com/Clinical-Genomics/genmod/pull/106))
- Released the constraint on Python 3.8 (collections, pkg_resources to importlib, tests) ([#142](https://github.com/Clinical-Genomics/genmod/pull/142))
- Update annotation examples ([#144](https://github.com/Clinical-Genomics/genmod/pull/144))
- Updated documentation with warning about compounds only being scored within the first family in the VCF ([#151](https://github.com/Clinical-Genomics/genmod/pull/151))
- Fixed sorting of variants ([#152](https://github.com/Clinical-Genomics/genmod/pull/152))
- genmod annotate for mitochondrial variants when using the `chrM` notation ([#157](https://github.com/Clinical-Genomics/genmod/pull/157))
- Use `uv` with `hatchling` to build and publish ([#104](https://github.com/Clinical-Genomics/genmod/issues/143))
- Fix linting issues ([#154](https://github.com/Clinical-Genomics/genmod/issues/154))

## [3.9]
- Fixed wrong models when chromosome X was named `chrX` and not `X` ([#135](https://github.com/Clinical-Genomics/genmod/pull/135))
- Added GitHub Actions workflows for automatic publishing to PyPI on release, and keep a changelog reminder ([#136](https://github.com/Clinical-Genomics/genmod/pull/136))
- Optional user defined threshold and penalty for compound scoring ([#138](https://github.com/Clinical-Genomics/genmod/pull/138))
- Update README with current github.io docs page ([#140](https://github.com/Clinical-Genomics/genmod/pull/140))

## [3.8.3]
- Fixed unstable compounds order in models output ([#134](https://github.com/Clinical-Genomics/genmod/pull/134))
- Added `six` to requirements.txt and setup.py ([#134](https://github.com/Clinical-Genomics/genmod/pull/134))

## [3.8.0]
- Rank score normalisation

## [3.7.4]

### Fixed
- OSError crashes due to socket address conflicts when using containers

## [x.x.x]

### Changed
- Adds support for escaped characters in FORMAT description header strings

### Fixed
- Documentation about cli options `strict` and `phased`
