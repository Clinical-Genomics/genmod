try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    

with open('README.txt') as file:
    long_description = file.read()

setup(name='genmod',
    version='0.6',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.commoonso/genmod',
    license = 'Modified BSD',
    install_requires=['ped_parser'],
    packages = ['genmod', 'genmod.utils', 'genmod.models', 'genmod.variants', 'genmod.vcf'],
    scripts = ['scripts/run_genmod.py'],
    keywords = ['inheritance', 'vcf', 'variants'],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)