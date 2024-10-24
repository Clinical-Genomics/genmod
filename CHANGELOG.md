# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

About changelog [here](https://keepachangelog.com/en/1.0.0/)

Please add a new candidate release at the top after changing the latest one. Feel free to copy paste from the "squash and commit" box that gets generated when creating PRs

Try to use the following format:

## [unreleased]
- Fixed wrong models when chromosome X was named `chrX` and not `X` 
- Added GitHub Actions workflows for automatic publishing to PyPI on release, and keep a changelog reminder ([#136](https://github.com/Clinical-Genomics/genmod/pull/136))
- Optional user defined threshold and penalty for compound scoring

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
