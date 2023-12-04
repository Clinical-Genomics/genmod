import sys
import os

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup
import pkg_resources

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
    version='3.8.2',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/genmod',
    license = 'MIT License',
    python_requires="~=3.8.0",
    install_requires=[
        'ped_parser >= 1.6.6',
        'pytabix >= 0.1',
        'pytest >= 7.3.1',
        'interval_tree >= 0.3.4',
        'click >= 8.1.3',
        'configobj >= 5.0.8',
        'intervaltree >= 3.1.0',
        'extract_vcf >= 0.5',
        'vcftoolbox >= 1.5.1'
    ],
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
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)
