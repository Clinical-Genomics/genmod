import sys

from codecs import (open, getreader)
import gzip
from multiprocessing import cpu_count
import click

variant_file = click.argument('variant_file', 
    type=click.Path(),
    metavar='<vcf_file> or -')

outfile = click.option('-o', '--outfile', 
    type=click.File('w'),
    help='Specify the path to a file where results should be stored.')

silent = click.option('-s', '--silent',
    is_flag=True,
    help='Do not print the variants.')

processes = click.option('-p', '--processes', 
    default=min(4, cpu_count()),
    help='Define how many processes that should be use for annotation.')

temp_dir = click.option('--temp_dir',
    type=click.Path(exists=True),
    help='Path to tempdir')

family_file = click.option('-f', '--family_file',
    type=click.File('r'),
    metavar='<ped_file>')

family_type = click.option('-t' ,'--family_type', 
    type=click.Choice(['ped', 'alt', 'cmms', 'mip']), 
    default='ped',
    help='If the analysis use one of the known setups, please specify which one.'
)


def get_file_handle(path):
    """Get a file handle"""
    if path == '-':
        if sys.version_info < (3,0):
            sys.stdin = getreader('utf-8')(sys.stdin)
        
        file_handle = sys.stdin
    
    elif path.endswith('.gz'):
        file_handle = gzip.open(path, 'r')
    
    else:
        file_handle = open(path, 'r')
    
    return file_handle