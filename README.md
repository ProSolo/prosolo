# ProSolo: bulk backing vocals for single cell solos

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/prosolo/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/prosolo/badges/downloads.svg)](http://bioconda.github.io/recipes/prosolo/README.html)

ProSolo is a variant caller for multiple displacement amplified DNA sequencing data from diploid single cells. It relies on a pair of samples, where one is from an MDA single cell and the other from a bulk sample of the same cell population, sequenced with any next-generation sequencing technology.

It uses an extension of the novel latent variable model of [Varlociraptor](https://github.com/varlociraptor/varlociraptor), that already integrates various levels of uncertainty. It adds a layer that accounts for amplification biases and errors of MDA, and thereby allows to properly asses the probability of having a variant or any other single cell event in the MDA single cell.

In the future, ProSolo will also implement indel calling, but currently, only the the `single-cell-bulk` subcommand with the `--omit-indels` flag is recommended.

## Installation

ProSolo is available via [Bioconda](https://bioconda.github.io), a distribution of bioinformatics software for the conda package manager.
[Bioconda can be set up](https://bioconda.github.io/user/install.html) in any Linux environment, even without admin rights.
With [Bioconda set up](https://bioconda.github.io/user/install.html), ProSolo can be installed into a dedicated conda environment via

    $ conda create -n prosolo prosolo


To then use it, you only need to activate its environment with

    $ conda activate prosolo

## Usage

For general usage, please issue `prosolo --help` on the command line. Help is also available for all subcommands listed there, e.g. `prosolo single-cell-bulk --help`.

### Variant Calling

To try out calling command syntax, please use the test data in the repo folder `tests/` as follows:
```
prosolo single-cell-bulk \
    --omit-indels \
    --candidates tests/candidates.bcf \
    --output test-out_omit-indels.bcf \
    tests/single-cell.bam \
    tests/bulk.bam \
    tests/ref.fa
```
As input, you always need:

1. Aligned reads of a single cell sample, in position sorted BAM or SAM format.
2. Aligned reads from a background bulk sample, in position sorted BAM or SAM format.
3. The reference FASTA file used in the alignment.
4. A file specifying which genomic sites to call and what possible alternative alleles to look at. This can be any VCF or BCF file, e.g. produced by `freebayes` or `samtools mpileup`.

If you want to parallelize variant calling, you can parallelize over genomic regions. Simply split your `--candidates` files into non-overlapping regions or even already parallelize the candidate file generation over genomic regions.

### Controlling the false discovery rate

To control the false discovery rate (FDR) for the minimal output of the above example, use its expected output data in the repo folder `tests/` as follows:
```
prosolo control-fdr \
    tests/expected-out_omit-indels.bcf \
    --events ADO_TO_REF,ADO_TO_ALT \
    --var SNV \
    --fdr 0.04
    --output ADO_fdr_0-04.bcf
```
In this case, we are jointly controlling the FDR for all `Events` that are allele dropouts in the single cell sample. For this set of `Events`, the above command will print a BCF file with all the entries with an estimated false discovery rate below `0.04`.

However, you can control the false discovery rate for any combination of events. The output of the `prosolo single-cell-bulk` command will always contain a `PROB_` field for every possible single cell event, and you can simply name the `Events` you jointly care about by omitting this prefix. E.g. to control the false discovery rate of the presence of an alternative allele, you would use `--events HET,HOM_ALT,ADO_TO_ALT,ADO_TO_REF,ERR_REF` to control for all the single cell events that imply the presence of an alternative allele.

### Using ProSolo in a pipeline

For systematic analyses of multiple single cell samples, we recommend using the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/latest/). To make setting up your own ProSolo workflow as easy as possible, we have created Snakemake wrappers for ProSolo's two subcommands, [`prosolo single-cell-bulk`](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/prosolo/single-cell-bulk.html) and [`prosolo control-fdr`](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/prosolo/control-fdr.html). For general usage of Snakemake wrappers, see the [Snakemake wrappers documentation](https://snakemake-wrappers.readthedocs.io/en/stable/index.html). And also have a look around for other wrappers: for most tools that you might need for preprocessing ProSolo input or for working with ProSolo results, there should already be Snakemake wrappers that are easy to integrate into a Snakemake pipeline.

## Authors

### Model

* [Johannes Köster](https://github.com/johanneskoester), [Louis Dijkstra](https://github.com/louisdijkstra) (see [varlociraptor](https://github.com/varlociraptor/varlociraptor))
* [David Lähnemann](https://github.com/dlaehnemann) (single cell whole genome amplification model, single cell & bulk joint calling model / event setup)
* [Alexander Schönhuth](https://github.com/aschoen) (latent variable model, single cell whole genome amplification model, single cell & bulk joint calling model / event setup)

### Implementation

* [Johannes Köster](https://github.com/johanneskoester) (see [varlociraptor](https://github.com/varlociraptor/varlociraptor))
* [David Lähnemann](https://github.com/dlaehnemann) (see [varlociraptor](https://github.com/varlociraptor/varlociraptor), all the implementation of the ProSolo CLI, originally based on the [prosic2 CLI](https://github.com/PROSIC/prosic2))

### Supervision

* Supervision of David Lähnemann: [Alice McHardy](https://github.com/alicemchardy) and [Alexander Schönhuth](https://github.com/aschoen)

### Affiliations

Affiliations during work on the project:

* [Life Sciences and Health group](https://www.cwi.nl/research/groups/life-sciences-and-health), Centrum Wiskunde & Informatica, Amsterdam, The Netherlands: Louis Dijkstra, Johannes Köster, Alexander Schönhuth
* [Computational Biology of Infection Research Group](https://www.helmholtz-hzi.de/en/research/research_topics/bacterial_and_viral_pathogens/computational_biology_of_infection_research/our_research/), Helmholtz Centre for Infection Research, Braunschweig, Germany: David Lähnemann, Alice McHardy
* [Algorithms for reproducible bioinformatics lab](https://koesterlab.github.io/), Institute of Human Genetics, University Hospital Essen, University of Duisburg-Essen, Germany: Johannes Köster, David Lähnemann
* [Department of Pediatric Oncology, Hematology, and Clinical Immunology](https://www.uniklinik-duesseldorf.de/en/unternehmen/kliniken/department-of-paediatric-oncology-haematology-and-immunology/), University Children’s Hospital, Medical Faculty, Heinrich Heine University, Düsseldorf, Germany: David Lähnemann
