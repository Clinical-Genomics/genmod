# Instructions for how to release a new version of Genmod

## Overview

Genmod uses `uv` and `hatchling` to build and publish the package to PyPI and Dockerhub.

## Genmod Git branching strategy

Genmod development is organised on a flexible Git ["Release Flow"][release_flow] branching system.
This more or less means that we make releases in release branches which corresponds to stable
versions of Genmod. The trunk branch is named `main`, and corresponds to the latest development version of Genmod.

The following steps are used to release a new version of Genmod.

1. Create a release branch with the release name, e.g. `version_3.12` and checkout

    ```bash
    git checkout -b version_3.12
    ```
   
2. Update version to, e.g. 3.12
   Either use the `uv version` tool to bump the version and update dependencies, or manually update the version in `genmod/pyproject.toml`.
   ```bash
      uv version --bump=minor
      git add genmod/pyproject.toml uv.lock
      git commit -m 'Bump version and upgrade dependencies'
   ```
   OR 
   - in `genmod/pyproject.toml`
     ```toml
     [project]
     name = "genmod"
     version = "3.12"
     ```
   ```bash
     git add genmod/pyproject.toml
     git commit -m 'Bump version'
   ```
   And update the `uv.lock` file to the latest compatible versions of dependencies
   This in particular also updates the genmod package verison. Note that it is now the same as in the previous step.
   ```bash
   uv lock --upgrade
   git add uv.lock
   git commit -m 'Upgrade dependencies'
   ```

4. Make sure CHANGELOG.md is up to date for the release. Add a new candidate release at the top after changing the latest one. 
 ```bash
   add CHANGELOG.md
    git commit -m 'Update CHANGELOG.md'
 ```
5. Push the release branch to github and create a pull request

    ```bash
    git push -u origin version_3.12
    ```
6. Thoroughly test any new features of Genmod, and ideally ensure they are covered by tests, manual example runs or external tests of failure points for bugfixes.
7. Solicit a review of the pull request and after getting it approved, merge it to main.
8. Draft a new release on GitHub, add some text - e.g. an abbreviated CHANGELOG - and release. 
9. The GitHub release action should take care of the rest (build and submit to PyPi and DockerHub, build and publish docs).
