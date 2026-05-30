#!/usr/bin/python3
"""
CLI entry point for Scaffold_builder.

Parses command-line arguments and delegates to the pipeline in
scaffold_builder.py. Configure logging here; keep core logic clean.
"""

import sys
import logging
import argparse

from utils import DEFAULTS
from utils import VERSION

from scaffold_builder import (build_parameters,
                               run)


def parse_args(argv=None):
    """Parse command-line arguments for Scaffold_builder.

    Args:
        argv (list[str] | None): Argument list; defaults to sys.argv[1:].

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="scaffold_builder",
        description=f"Scaffold_builder v{VERSION}: Combining de novo and reference-guided assembly.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-q", required=True, metavar="QUERY",
                        help="Query contigs in FASTA format.")
    parser.add_argument("-r", required=True, metavar="REFERENCE",
                        help="Reference genome in FASTA format.")
    parser.add_argument("-p", default=DEFAULTS["-p"], metavar="PREFIX",
                        help=f"Output file prefix (default: {DEFAULTS['-p']}).")
    parser.add_argument("-t", type=int, default=DEFAULTS["-t"], metavar="INT",
                        help=f"Terminus length to align (default: {DEFAULTS['-t']} nt).")
    parser.add_argument("-i", type=int, default=DEFAULTS["-i"], metavar="INT",
                        help=f"Minimum identity for merging (default: {DEFAULTS['-i']}%%).")
    parser.add_argument("-a", type=int, default=DEFAULTS["-a"], metavar="INT",
                        help=f"Ambiguous mapping threshold (default: {DEFAULTS['-a']}%%).")
    parser.add_argument("-b", type=int, default=DEFAULTS["-b"], choices=[0, 1],
                        help=f"Rearrangement behaviour: 0=end-to-end, 1=new scaffold (default: {DEFAULTS['-b']}).")
    parser.add_argument("-g", type=int, default=DEFAULTS["-g"], metavar="INT",
                        help=f"Maximum gap length (default: {DEFAULTS['-g']} nt).")
    return parser.parse_args(argv)


def main():
    """Configure logging, parse arguments, and run the pipeline."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    parameters = build_parameters({
        "-q": args.q, "-r": args.r, "-p": args.p,
        "-t": args.t, "-i": args.i, "-a": args.a,
        "-b": args.b, "-g": args.g,
    })
    try:
        run(parameters)
        print("Scaffolded :)")
    except Exception as exc:
        logging.error("%s", exc)
        sys.exit(1)
