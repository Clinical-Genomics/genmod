try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import pkg_resources

# For making things look nice on pypi:
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except (IOError, ImportError, RuntimeError):
    long_description = 'Tool for annotating patterns of genetic inheritance in Variant Call Format (VCF) files.'


# with open('README.txt') as file:
#     long_description = file.read()

setup(name='genmod',
    version='2.0.4',
    description='Annotate genetic inheritance models in variant files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/genmod',
    license = 'MIT License',
    install_requires=[
        'ped_parser >= 1.2.3',
        'vcf_parser == 1.1.4',
        'pytabix',
        'pytest',
        'interval_tree',
        'click',
        'configparser',
        'logbook'
    ],
    packages = [
        'genmod',
        'genmod/commands',
        'genmod/utils',
        'genmod/errors',
        'genmod/models',
    ],
    package_data = {
        'genmod': ['annotations/*.db']
    },
    scripts = [
        'scripts/genmod'
    ],
    # entry_points= { "console_scripts" : [
    #     "run_genmod = scripts.run_genmod:main",
    # ]},
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