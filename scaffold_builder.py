#!/usr/bin/python3
"""
Scaffold_builder pipeline orchestration.

Provides build_parameters(), run_nucmer(), and run() which wire together
the pure utility functions from utils.py into the full scaffolding pipeline.
"""

import os
import logging
import subprocess

from shutil import which
from pathlib import Path

from utils import (DEFAULTS,
                   VERSION,
                   clean_coords,
                   coord2hash,
                   fasta2hash,
                   sort_mapping)


def build_parameters(user_args):
    """Merge user-supplied arguments into the default parameters dict.

    Args:
        user_args (dict): Flag → value overrides (values may be str or int).

    Returns:
        dict: Full parameters dict with runtime state initialised.
    """
    params = {
        **DEFAULTS,
        "begin_end": 0,
        "ambiguous": [],
        "subSet": [],
        "rawMapping": "",
        "gaps": [],
        "overlapOver": 0,
        "overlapBelow": 0,
    }
    params.update(user_args)
    for flag in ("-t", "-g", "-i", "-a", "-b"):
        params[flag] = int(params[flag])
    return params


def run_nucmer(parameters):
    """Run nucmer and show-coords to produce the .coords alignment file.

    Args:
        parameters (dict): Requires '-r', '-q', '-p' keys.

    Raises:
        FileNotFoundError: If nucmer or show-coords are absent from PATH.
        subprocess.CalledProcessError: If either tool exits non-zero.
    """
    nucmer_bin = which("nucmer")
    show_coords_bin = which("show-coords")
    if not nucmer_bin:
        raise FileNotFoundError("nucmer not found on PATH")
    if not show_coords_bin:
        raise FileNotFoundError("show-coords not found on PATH")

    prefix = parameters["-p"]
    Path(prefix + "_overlap_alignment").mkdir(exist_ok=True)

    subprocess.run(
        [nucmer_bin, "-b", "1500", "-g", "500",
         parameters["-r"], parameters["-q"], "-p", prefix],
        check=True,
    )
    with open(prefix + ".coords", "w") as fh:
        subprocess.run([show_coords_bin, prefix + ".delta"], stdout=fh, check=True)
    os.remove(prefix + ".delta")
    logging.info("nucmer alignment complete: %s.coords", prefix)


def run(parameters):
    """Execute the full Scaffold_builder pipeline.

    Runs nucmer, parses coordinates, removes ambiguous mappings, and
    scaffolds the query contigs against the reference.

    Args:
        parameters (dict): Fully populated dict from build_parameters().
    """
    run_nucmer(parameters)
    hash_sequences = fasta2hash(parameters["-q"])
    coords = coord2hash(parameters["-p"] + ".coords", parameters, hash_sequences)
    cleaned = clean_coords(coords, parameters, hash_sequences)
    sort_mapping(cleaned, parameters, hash_sequences)
    logging.info("Scaffolded successfully.")
