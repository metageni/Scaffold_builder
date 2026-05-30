![](logo/scaffold_builder_logo.png "Logo")

#### Scaffold_builder: Combining de novo and reference-guided assembly with Scaffold_builder
* [Installation](#installation)
    * [dependencies](#dependencies)
    * [pip](#pip)
    * [bioconda](#bioconda)
    * [git](#git)
* [Usage](#usage)
* [Testing](#testing)
* [Citing](#citing)

> **v3.0.0** — Complete Python 3 rewrite. Produces byte-for-byte identical scaffold
> output to v2.2. Fixes N50 calculation and average-length rounding bugs present in
> the original. See [CHANGELOG.md](CHANGELOG.md) for details.

## Installation
### Dependencies
- [Python 3.8+](https://www.python.org/downloads/)
- [click](https://click.palletsprojects.com/)
- [MUMmer4 (nucmer)](https://github.com/mummer4/mummer)


### Bioconda
Install Scaffold_builder and all dependencies (including MUMmer4) via [conda](https://conda.io):

    conda create -n scaffold_builder -c bioconda -c conda-forge scaffold-builder
    conda activate scaffold_builder

### Git

    git clone https://github.com/metageni/Scaffold_builder.git
    cd Scaffold_builder
    pip install .
    scaffold_builder -q [QUERY] -r [REFERENCE]

## Usage

    scaffold_builder -q query_contigs.fna -r reference_genome.fna -p output_prefix [-t] [-i] [-a] [-b] [-g]

    -q  fasta file of contigs
        Required. Query contigs in FASTA format.

    -r  fasta file containing reference genome
        Required. Reference genome in FASTA format.

    -p  prefix for output files  [default: Scaffold]

    -t  length of terminus that will be aligned  [default: 300 nt]
        Scaffold_builder checks whether the termini of adjacent contigs are
        homologous by aligning them with Needleman-Wunsch and merges them if so.

    -i  minimum identity for merging contigs  [default: 80%]
        Non-identical aligned nucleotides are represented with IUPAC ambiguity codes.

    -a  ambiguous mapping threshold  [default: 95%]
        Contigs mapping to more than one location over this percentage of their
        length are excluded from scaffolding.

    -b  rearrangement behaviour  [default: 0]
        0: place end-to-end
        1: create new scaffold sequence

    -g  maximum gap length allowed  [default: 5000 nt]
        Gaps larger than this value split the scaffold into a new sequence.

## Citing
Scaffold_Builder was written by Genivaldo G. Z. Silva. Feel free to [contact me](mailto:genivaldo.gueiros@gmail.com)

If you use Scaffold_Builder, please cite it:

	Silva GG, Dutilh BE, Matthews TD, Elkins K, Schmieder R, Dinsdale EA, Edwards RA.
	Combining de novo and reference-guided assembly with Scaffold_builder,
	Source Code for Biology and Medicine 2013.
