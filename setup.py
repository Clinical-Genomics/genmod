try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    
# For making things look nice on pypi:
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError):
    long_description = 'Tool for annotating patterns of genetic inheritance in Variant Call Format (VCF) files.'


# with open('README.txt') as file:
#     long_description = file.read()

setup(name='genmod',
    version='1.2.2',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/genmod',
    license = 'MIT License',
    install_requires=['ped_parser', 'vcf_parser', 'pysam', 'pytest', 'interval_tree'],
    packages = ['genmod', 'genmod.utils', 'genmod.models', 'genmod.variants'],
    scripts = ['scripts/run_genmod.py'],
    keywords = ['inheritance', 'vcf', 'variants'],
    classifiers = [
        "Programming Language :: Python",
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