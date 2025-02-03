#!/usr/bin/env python
# encoding: utf-8
"""
run_genmod.py

Script for annotating genetic models in variant files.

Created by Måns Magnusson on 2014-01-21.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function

import click

# , sort, annotate, analyze, summarize_variants, score_variants)
from genmod import __version__

from . import (
    annotate_variant_command,
    filter_command,
    models_command,
    score_command,
    score_compounds_command,
    sort_command,
)


def print_version(ctx, param, value):
    """Callback function for printing version and exiting
    Args:
        ctx (object) : Current context
        param (object) : Click parameter(s)
        value (boolean) : Click parameter was supplied or not
    Returns:
        None:
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo("genmod version: " + __version__)
    ctx.exit()


###         This is the main script         ###


@click.group()
@click.option("--version", is_flag=True, callback=print_version, expose_value=False, is_eager=True)
@click.option(
    "-l",
    "--logfile",
    type=click.Path(exists=False),
    help="Path to log file. If none logging is printed to stderr.",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=0,
    help="Increase output verbosity. Can be used multiple times, eg. -vv",
)
@click.pass_context
def cli(context, logfile, verbose):
    """Tool for annotating and analyzing genetic variants in the vcf format.\n
    For more information, please run:
    genmod COMMAND --help \n
    """
    from genmod import logger
    from genmod.log import LEVELS, init_log

    loglevel = LEVELS.get(min(verbose, 2), "WARNING")

    init_log(logger, logfile, loglevel)


cli.add_command(sort_command)
cli.add_command(models_command)
cli.add_command(score_command)
cli.add_command(score_compounds_command)
cli.add_command(annotate_variant_command)
cli.add_command(filter_command)


if __name__ == "__main__":
    cli()
