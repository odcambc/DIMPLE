[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/coywil26/DIMPLE/blob/master/DIMPLE.ipynb)

# DIMPLE: Deep Indel Missense Programmable Library Engineering

## Protein domain insertion via programmed oligo libraries

A Python script for generating oligo libraries and PCR primers for Deep Mutational Scanning library generation incorporating indel variation.

Take a look at the protocol on [protocols.io](https://www.protocols.io/view/dimple-library-generation-and-assembly-protocol-rm7vzy7k8lx1) for more information on generating and assembling libraries as well.

Our pipeline for running and analyzing screens with DIMPLE libraries, [dumpling](https://github.com/odcambc/dumpling), is available as well.

Note: This is the active repository for DIMPLE development. The archived repository containing the code used in the publication is [here](https://github.com/odcambc/DIMPLE), and is also archived at [Zenodo](https://zenodo.org/records/7574261).

# Quick start

The fastest way to run DIMPLE is a [Google Colab](https://colab.research.google.com/github/coywil26/DIMPLE/blob/master/DIMPLE.ipynb) notebook. Follow the prompts and generate a library.

You can also run the DIMPLE tool as a web app at [Hugging face](https://huggingface.co/spaces/cbcbcbcbcb/dimple-demo).

# Installation

DIMPLE includes a command-line version and a GUI version. Both of them come with an install.

## Using uv (recommended)

DIMPLE is an application you download and run, not a library you install. `uv`
reads `pyproject.toml` + `uv.lock`, creates a local `.venv` with the locked
dependencies, and runs DIMPLE.

From a fresh clone, a single command bootstraps everything and runs:

```{bash}
uv run python run_dimple.py -h
# or run the gui:
uv run python run_dimple_gui.py
```

## Without uv

You can use conda:

```{bash}
conda create -n dimple python=3.12
conda activate dimple
pip install biopython==1.84 numpy==1.26.4 pydna==5.4.0
python run_dimple.py -h
# or run the gui:
python run_dimple_gui.py
```

Or you can manually install the dependencies and then run DIMPLE from the command line like above.

# Inputs

## Target gene file

Targeted genes should be supplied in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). To allow DIMPLE to check for nonspecific amplification, include the entire plasmid sequence of the library generation construct in the file.

The ORF can be specified in the fasta header for each target gene. If desired, the header should include the start and end positions of the gene in the plasmid, as follows:

```{text}
>gene1 start:10 end:100
ATGTT...
```

The start position should be the first base of the first codon, and the end position should be the last base of the last codon. Otherwise specify the ORF in the command line.

# Running DIMPLE

## Colab version

Using the [Google Colab](https://colab.research.google.com/github/coywil26/DIMPLE/blob/master/DIMPLE.ipynb) notebook, follow the prompts and explanations. Also check the options below for additional usage.

## Local version

We have supplied two methods to run DIMPLE: a command-line version, and a GUI.
Both offer the same functionality, but the GUI is more user-friendly.

### GUI usage

Start the GUI with the following command:

```{bash}
uv run python run_dimple_gui.py
# or, without uv:
python run_dimple_gui.py
```

![DIMPLE_GUI](DIMPLE/data/DIMPLE_GUI.png)

The following are required:

- A directory for results (working directory)
- Target gene file (see below for format requirements)
- One or more of the mutations to make to the target gene

Supply options, then generate library by pressing 'Run DIMPLE' button.

### Command-line usage

See a description of options for command-line version:

```{bash}
uv run python run_dimple.py -h
# or, without uv:
python run_dimple.py -h
```

Full list of options:

<details open>
<summary> Full list of command-line options</summary>

```{text}
options:
  -h, --help            show this help message and exit
  -wDir WDIR            Working directory for fasta files and output folder
  -geneFile GENEFILE    Input all gene sequences including backbone in a fasta format. Place all in one fasta file.
                        Name description can include start and end points (>gene1 start:1 end:2)
  -handle HANDLE        Genetic handle (linker) sequence that -dis inserts at every position. Important for defining
                        the linker. Currently uses BsaI (4 base overhang), but this can be swapped for SapI (3 base
                        overhang).
  -dis                  Domain-insertion scan: insert the -handle sequence at every position in the ORF.
  -matchSequences       Find similar sequences between genes to avoid printing the same oligos multiple times.
                        Default: No matching
  -oligoLen OLIGOLEN    Total synthesized oligo length (nt), i.e. the full oligo your synthesis vendor makes; the
                        usable mutagenic fragment is this minus fixed primer/cutsite overhead
  -fragmentLen FRAGMENTLEN
                        Maximum length of gene fragment
  -overlap OVERLAP      Enter number of bases to extend each fragment for overlap. This will help with insertions
                        close to fragment boundary
  -DMS, -include_substitutions
                        Run a deep mutational (substitution) scan. -include_substitutions is an alias, matching the
                        GUI's checkbox name for the same toggle.
  -custom_mutations CUSTOM_MUTATIONS
                        Path to file that includes custom mutations with the format position:AA
  -usage USAGE          Default is "human" or "ecoli", or pass a path to a codon usage file
  -insertions INSERTIONS [INSERTIONS ...]
                        Enter a list of insertions (nucleotides) to make at every position. Note, you should enter
                        multiples of 3 nucleotides to maintain reading frame
  -deletions DELETIONS [DELETIONS ...]
                        Enter a list of deletions (number of nucleotides) to symmetrically delete (it will make
                        deletions in multiples of 2x). Note you should enter multiples of 3 to maintain reading frame
  -barcode_start BARCODE_START
                        To run DIMPLE multiple times, you will need to avoid using the same barcodes. This allows you
                        to start at a different barcode.
  -restriction_sequence RESTRICTION_SEQUENCE
                        Recommended using BsmBI - CGTCTC(G)1/5 or BsaI - GGTCTC(G)1/5. Do not use N
  -avoid_sequence AVOID_SEQUENCE [AVOID_SEQUENCE ...]
                        Avoid these sequences in the backbone - BsaI and BsmBI. For multiple sequences use a space
                        between inputs. Example -avoid_sequence CGTCTC GGTCTC
  -include_stop_codons  Include stop codons in the list of scanning mutations.
  -include_synonymous   Include synonymous codons in the list of scanning mutations.
  -make_double          Make each combination of mutations within a fragment
  -maximize_nucleotide_change
                        Maximize the number of nucleotide changes in each codon for easier detection in NGS
  -seed SEED            Seed for random number generation (default: 1848 for reproducible output; pass any integer,
                        including 0, to vary the library)
  --non_interactive     Run without interactive prompts (fail fast on ambiguous ORF selection).
  --orf_index ORF_INDEX
                        Preferred ORF index for non-interactive ORF selection.
  --link_policy {prompt,always,never}
                        Policy for linking matched genes during align_genevariation.
  --breaksite_change_policy {prompt,warn,error}
                        Policy when breaksite endpoints change outside DMS mode.
```

</details>

# DIMPLE config

DIMPLE is flexible enough to create libraries with different types of mutations.

- Stop codons
- Synonymous codons
- Non-synonymous codons
- Insertions
- Deletions
- Domain insertion handles
- Double mutations (in subpools)

There are a few important options.

- Oligo length: this is the desired length of the **synthesized oligos**
- Type IIS restriction sequence: this is the recognition sequence of the IIS enzyme that will be used for assembly
- Maximize nucleotide change: this picks codons with the maximum number of base changes to potentially improve variant detection
- Double fragment per oligo: this tries to pack two fragments into a single oligo for ordering efficiency. It is currently not implemented.

# DIMPLE output

## Primer lists

- List of all primers: All_Primers.fasta
- List of primers for amplifying backbone: Example_DMS_Gene_Primers.fasta
- List of primers for amplifying oligos: Example_DMS_Oligo_Primers.fasta

## Oligo lists

- List of all oligos in the pool: All_Oligos.fasta
- List of gene-specific oligos: Example_DMS_Oligos.fasta

## Designed variants

This is a list of the designed variants that can be used directly as input to our experimental pipeline [dumpling](github.com/odcambc/dumpling)

## Examples

Example output files are located in the `examples` directory.

# Running tests

To test DIMPLE, run the following command from the root directory:

  ```{bash}
uv run pytest
# or, without uv:
pytest
```

This should pass without any errors. If you encounter any issues, please open an issue on the GitHub repository.

# Troubleshooting and known issues

## Version incompatibilities

DIMPLE has been tested on Python version 3.12. Biopython is currently incompatible with Python 3.13 in some cases, and we recommend using Python 3.12 for now.

## ORF issues with non-interactive runs

DIMPLE auto-detects ORFs by scanning all six reading frames and looking for stretches longer than 100 amino acids. In an interactive run this prompts the user to pick one if the ORF isn't specified in the header. Colab runs may not prompt, so you should set the ORF in the FASTA header:

```{text}
>gene1 start:10 end:100
ATGTT...
```

If no ORF is provided and DIMPLE has an issue detecting the ORF, it will raise a `ValueError`:

- *No ORF candidates found* → no ≥100-aa frame matched. Add explicit `start:`/`end:` header.
- *Multiple ORF candidates found* → multiple ≥100-aa frames matched. Add `start:`/`end:` header or `--orf_index N`.
- *Preferred ORF index N is out of range* → `--orf_index` value larger than the detected count.

The Colab notebook above defaults to `non_interactive=True` and has an `orf_index` form field for this reason.

## Failure to generate a library

Sometimes DIMPLE won't be able to find a set of oligos and primers that satisfy the settings. This is more likely to happen with short oligos and large indels. It can also happen if double mutations are selected.

A `PCR not specific` error usually indicates an indel problem, but a `Primers no longer bind` issue is a double mutation problem. Try adjusting mutation settings or oligo length.

## GUI seems to hang or exits with EOFError

If you don't supply an ORF or if it's misformatted, the GUI will ask for help identifying it on the console: if you launched the DIMPLE GUI from the command line, it will be there. Try making sure your ORF is correct in the input fasta.

## Frameshifts and indel problems

Insertions and deletions are specified in **nucleotides**: codon-wise deletions would be deletions of 3,6,9,... and **not** 1,2,3...

DIMPLE can create frameshift libraries, but the variants won't be translated in the output files.

## DIMPLE runs out of barcodes

DIMPLE uses a set of unique barcodes to orthogonally amplify subpools. Sometimes a library will end up needing more than the 3000 it has. You might have to change the number of mutations or genes in the input.

## Barcode clashes

By default, DIMPLE starts generating a library from the first barcode. That means that libraries generated from two runs will have the **same** barcodes by default - and the primers won't be specific as a result. If you want to create an oligo pool for a number of different genes, make sure you supply them as a single combined fasta. Setting the starting barcode is possible but not recommended.

## Double fragments per oligo (doublefrag) is disabled

 This tries to pack two fragments into a single oligo for ordering efficiency and would only be recommended for rare cases with extra testing. It is currently not implemented.

# Citing DIMPLE

If you found DIMPLE useful, feel free to cite the publication describing it:

- Preprint: [Macdonald et al., 2022](https://doi.org/10.1101/2022.07.26.501589)
- Published: [Macdonald et al., 2023](https://doi.org/10.1186/s13059-023-02880-6)

# License

This code is licensed under the terms of the MIT license: [License](LICENSE)

# Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.waymentsteelelab.org).
