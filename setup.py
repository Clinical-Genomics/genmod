import io
import sys
import os

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

here = os.path.join(os.path.abspath(os.path.dirname(__file__)))

def parse_reqs(req_path = "./requirements.txt"):
    """Recursively parse requirements from nested pip files."""
    install_requires = []
    with io.open(os.path.join(here, req_path), encoding="utf-8") as handle:
        # remove comments and empty lines
        lines = (line.strip() for line in handle if line.strip() and not line.startswith("#"))

        for line in lines:
            # check for nested requirements files
            if line.startswith("-r"):
                # recursively call this function
                install_requires += parse_reqs(req_path=line[3:])
            else:
                # add the line as a new requirement
                install_requires.append(line)

    return install_requires

about = {}
with open(os.path.join(here, "genmod/__version__.py")) as f:
    exec(f.read(), about)

# What packages are required for this module to be executed?
REQUIRED = parse_reqs()

# Shortcut for building/publishing to Pypi
if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist bdist_wheel upload')
    sys.exit()

# For making things look nice on pypi:
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError, RuntimeError):
    long_description = 'Tool for annotating patterns of genetic inheritance in Variant Call Format (VCF) files.'

setup(name='genmod',
    version='3.9',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/Clinical-Genomics/genmod',

    license = 'MIT License',
    install_requires=REQUIRED,
    packages=find_packages(
        exclude=('tests*', 'docs', 'examples', 'configs')
    ),
    # package_data = {
    #     'genmod': ['annotations/*']
    # },
    include_package_data = True,


    entry_points= { "console_scripts" : [
        "genmod = genmod.commands.base:cli",
        ]
    },

    keywords = ['inheritance', 'vcf', 'variants'],
    classifiers = [
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)
