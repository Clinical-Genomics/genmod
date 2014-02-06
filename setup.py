try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    

with open('README.txt') as file:
    long_description = file.read()

setup(name='genmod',
      version='0.2',
      description='Annotate genetic inheritance models in variant files',
      author = 'Mans Magnusson',
      author_email = 'mans.magnusson@scilifelab.se',
      url = 'http://github.commoonso/genmod',
      license = 'Modified BSD',
      packages = ['genmod', 'genmod.utils', 'genmod.models', 'genmod.variants', 'genmod.family'],
      scripts = ['scripts/run_genmod.py'],
      long_description = long_description,
)