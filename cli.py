#!/usr/bin/python3
"""
CLI entry point for Scaffold_builder.

Defines the command-line interface using click and delegates to the
scaffolding pipeline in scaffold_builder.py. Logging is configured here.
"""

import sys
import logging

import click

from utils import DEFAULTS
from utils import VERSION

from scaffold_builder import (build_parameters,
                               run)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(VERSION, "-v", "--version")
@click.option("-q", required=True, type=click.Path(exists=True),
              help="Query contigs in FASTA format.")
@click.option("-r", required=True, type=click.Path(exists=True),
              help="Reference genome in FASTA format.")
@click.option("-p", default=DEFAULTS["-p"], show_default=True,
              help="Output file prefix.")
@click.option("-t", default=DEFAULTS["-t"], show_default=True, type=int,
              help="Terminus length to align (nt).")
@click.option("-i", default=DEFAULTS["-i"], show_default=True, type=int,
              help="Minimum identity for merging overlaps (%%).")
@click.option("-a", default=DEFAULTS["-a"], show_default=True, type=int,
              help="Ambiguous mapping threshold (%%).")
@click.option("-b", default=str(DEFAULTS["-b"]), show_default=True,
              type=click.Choice(["0", "1"]),
              help="Rearrangement behaviour: 0=end-to-end, 1=new scaffold.")
@click.option("-g", default=DEFAULTS["-g"], show_default=True, type=int,
              help="Maximum gap length before scaffold break (nt).")
def main(q, r, p, t, i, a, b, g):
    """Scaffold_builder: combine de novo and reference-guided assembly.

    Scaffolds query contigs (-q) against a reference genome (-r) using
    MUMmer (nucmer) alignments.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parameters = build_parameters({
        "-q": q, "-r": r, "-p": p,
        "-t": t, "-i": i, "-a": a,
        "-b": int(b), "-g": g,
    })
    try:
        run(parameters)
        click.echo("Scaffolded :)")
    except Exception as exc:
        logging.error("%s", exc)
        sys.exit(1)

if __name__ == "__main__":
    main()
