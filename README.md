# Processing potential variant calls in single cells (ProSolo)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/prosolo/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/prosolo/badges/downloads.svg)](http://bioconda.github.io/recipes/prosolo/README.html)

ProSolo is a variant caller for single cell data from whole genome amplification with multiple displacement amplification (MDA). It relies on a pair of samples, where one is from an MDA single cell and the other from a bulk sample of the same cell population, sequenced with any next-generation sequencing technology.

It uses an extension of the novel latent variable model of [libprosic](https://github.com/PROSIC/libprosic), that already integrates various levels of uncertainty. It adds a layer that accounts for amplification biases (and errors) of MDA, and thereby allows to properly asses the probability of having a variant in the MDA single cell.

In the future, ProSolo will also implement indel calling, but currently, only the the `single-cell-bulk` subcommand with the `--omit-indels` flag is recommended.

## Installation

ProSolo is available via [Bioconda](https://bioconda.github.io), a distribution of bioinformatics software for the conda package manager.
[Bioconda can be set up](https://bioconda.github.io/#using-bioconda) in any Linux environment, even without admin rights.
With [Bioconda set up](https://bioconda.github.io/#using-bioconda), ProSolo can be installed via

	$ conda install prosolo

## Usage

### Variant Calling

To try out calling command syntax, please use the test data in the repo folder `tests/` as follows:
```
prosolo single-cell-bulk \
		--omit-indels \
    --candidates tests/candidates.bcf \
    --output test-out_omit-indels.bcf \
    --sc-isize-mean 12 \
    --sc-isize-sd 1 \
    tests/single-cell.bam \
    tests/bulk.bam \
    tests/ref.fa
```

### Controlling the false discovery rate

To control the false discovery rate (FDR) for the minimal output of the above example, use its expected output data in the repo folder `tests/` as follows:
```
prosolo control-fdr \
    tests/expected-out_omit-indels.bcf \
    --events ADO_TO_REF,ADO_TO_ALT \
    --var SNV \
    --method ev
```
In this case, we are jointly controlling the FDR for all `Events` that are allele dropouts in the single cell sample. For this set of `Events`, the above command will print probability cutoffs for different false discovery rates, e.g. `0.1819080263` for a false discovery rate of `0.04`. You can then use this to filter variants on the joint probabilities of the respective events:
```
prosolo apply-fdr \
    tests/expected-out_omit-indels.bcf \
    --events ADO_TO_REF,ADO_TO_ALT \
    --var SNV \
    --threshold 0.1819080263 \
		--output ADO_fdr_0-04.bcf
```

# Authors

* Original model: [Louis Dijkstra](https://github.com/louisdijkstra)
* Extended model and implementation (libprosic and PROSIC2): [Johannes Köster](https://johanneskoester.bitbucket.org)
* Extended model and extended implementation (libprosic and ProSolo): [David Laehnemann]()
