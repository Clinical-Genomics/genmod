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
    version='3.5.3',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/genmod',
    license = 'MIT License',
    install_requires=[
        'ped_parser >= 1.6.2',
        'pytabix',
        'pytest',
        'interval_tree >= 0.3.2',
        'click',
        'configobj',
        'intervaltree',
        'extract_vcf >= 0.4.2'
    ],
    packages=find_packages(
        exclude=('tests*', 'docs', 'examples', 'configs')
    ),
    package_data = {
        'genmod': ['annotations/*.db']
    },
    entry_points= { "console_scripts" : [
        "genmod = genmod.commands.run_genmod:cli",
        ]
    },
    keywords = ['inheritance', 'vcf', 'variants'],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)