# Software Manual for SatuTe

## 1. Introduction

### Overview

Welcome to the SatuTe manual. This document provides comprehensive information on how to install, configure, and use SatuTe effectively.

SatuTe (**Satu**ration **Te**st) is a Python-based tool designed to test for phylogenetic information in phylogenetic analyses. The absence of phylogenetic information can be considered saturation. For two sequences, saturation occurs when multiple substitutions obscure true genetic distances, potentially leading to artifacts and errors. SatuTe provides a new measure that generalizes the concept of saturation between two sequences to a theory of saturation between subtrees. The implemented test quantifies whether the given alignment provides enough phylogenetic information shared between two subtrees connected by a branch in a phylogeny.

This enables the detection of branch saturation and assesses the reliability of inferred phylogenetic trees and the data from which they are derived in phylogenetic reconstruction.

### Requirements

The minimal input of SatuTe is a multiple sequence alignment, a model of sequence evolution with its parameters, and a phylogenetic tree. SatuTe parses these essential pieces of information from the output of [IQ-Tree](http://www.iqtree.org/), an efficient software for phylogenomic inference. While we strongly recommend running IQ-Tree separately with customized options, SatuTe can also use an IQ-Tree executable to generate any missing information using default settings.

**Technical Requirements:**

- Python: 3.6 or higher
- IQ-Tree: 2.2.2.3 or higher

### Main Workflow

The main function of SatuTe operates as follows: Given the required input, SatuTe first calculates the spectral decomposition of the rate matrix and determines the likelihood vectors for each node in the tree. It then performs the test for phylogenetic information on a user-selected branch or on each branch of the tree, as described in the relevant literature. The program outputs the test results and its components in different CSV files and a Nexus file (see [SatuTe Manual PDF](./docs/SatuTe_Manual.pdf)).

![Theoretical Foundation of SatuTe](./docs/figure_1_2024_08_02.png)
**Figure:** **Theoretical Foundation of SatuTe**

**a**, In an alignment of the sequences from five taxa, each column represents a pattern $\partial$.

**b**, The branch $AB$ splits the five-taxon tree into subtrees $\mathbb{T}_A$ with sequences $S1,S3,S4$ and subtree $\mathbb{T}_B$ with sequences $S2, S5$. Additionally, branch $AB$ splits each pattern, such as $\partial = \text{CTTCT}$, into subpatterns ${\partial A=\text{CTC}}$ and ${\partial B=\text{TT}}$. These subpatterns are used to compute the likelihood vectors ${{L}(\partial A)}$ and ${{L}(\partial B)}$, which are then used to compute the coefficient $C^\partial_1$.

**c**, The average coherence $\hat{C}_1$ is approximately normal distributed. Under the null hypothesis $H_0$ that the alignment was generated independently by subtrees $\mathbb{T}_A$ and $\mathbb{T}_B$, the expectation is zero, while its variance $\hat{\sigma}^2_1/n$ is easier to estimate. Combining this, we can compute the SatuTe z-score $Z$.

**d**, If $Z>z_\alpha$, we reject the null hypothesis $H_0$ at the significance level $\alpha$ and say that the alignment is phylogenetic informative for this branch. Otherwise, the alignment is saturated, meaning that it does not reflect the existence of branch $AB$.

In cases where a model of rate heterogeneity is used, SatuTe assigns each site to the rate category with the highest posterior probability. The alignment is then split by SatuTe. For each category, the test for branch saturation is employed on the rescaled phylogenetic tree given the subalignment.

## 2. Installation

Satute is available as a python package from pypi and can be normally installed via pip.
We recommend to use [pipx](https://pipx.pypa.io/stable/) to install Satute as a standalone command line tool. Using pipx ensures that Satute and its dependencies are installed in an isolated environment, minimizing potential conflicts with other Python packages on your system.

### Prerequisites

- Python 3.6 or higher
- `pipx` (Python package installer)

### Install Satute using pipx

1. **Install pipx:**  If you don't have pipx installed, you can install it using pip:

    ```bash
    pip install pipx
    ```

    After installation, ensure pipx is set up correctly:

    ```bash
    pipx ensurepath
    ```

2. **Install Satute using pipx:**  Once pipx is installed, you can use it to install Satute:

    ```bash
    pipx install satute
    ```

For more detailed instructions and information about pipx, refer to the [official pipx documentation](https://pipxproject.github.io/pipx/).

### Verifying the Installation

After the installation is complete, you can verify that Satute has been installed correctly by running the following command:

```bash
satute version
```

You should see the version number of Satute printed on the screen, confirming that the installation was successful.

## 3. Getting Started

If you have a previous output from an [IQ-Tree](http://www.iqtree.org/) run, such as one from the Webserver, you can specify the directory using the `-dir` option for an easy start.

```bash
satute -dir ./examples/example_dna/dir_ML_tree
```

or as alternative: Given a path to an IQ-Tree executable, SatuTe runs IQ-Tree with default options to generate the required data for the test of saturation, namely a multiple sequence alignment (MSA), a model of sequence evolution with its parameters, and a phylogenetic tree.

```bash
satute  -msa ./examples/data/example_dna.fasta -iqtree path_to_iqtree_exe
```
  
In our examples folder, we provide files for showcasing SatuTe on amino acid alignments in the subdirectory ‘./examples/examples_aa’. For more detailed information and three explained examples, please refer to the [SatuTe Manual PDF](./docs/SatuTe_Manual.pdf).

## 4. Options of SatuTe

| Option | Description  |
| ------------ | ----------------------------------------------------------- |
| `-dir <directory_path>` | Path to the input directory containing IQ-Tree output files. Use this option when you've already run IQ-Tree and want to avoid rerunning it. The directory should contain essential IQ-Tree output files including the .iqtree file, tree file(s), and possibly a .siteprob file. |
| `-alpha <significance_level>` | Significance level for the saturation test. The default value is 0.05, indicating a 5% significance level. |
| `-msa <msa_file_path>` | Path to the Multiple Sequence Alignment (MSA) file you wish to analyze. The MSA can be in FASTA, NEXUS, PHYLIP, or TXT format. |
| `-model <evolution_model>` | Indicates the model of sequence evolution. Common models include GTR, HKY, etc. You can also specify rate heterogeneity and other model extensions, like +G4 for gamma-distributed rates. |
| `-tree <tree_file_path>` | Path to the input tree file in Newick or Nexus format. This tree will be used as the basis for the saturation analysis. |
| `-iqtree <iqtree_path>` | Specifies the path to the IQ-Tree executable. If IQ-Tree is installed system-wide, just providing the executable name (iqtree or iqtree2) will suffice. |
| `-add_iqtree_options <additional_option>` | Specify additional options for the IQ-Tree run, if necessary. |
| `-edge <edge_name>` | Specify a branch or edge name to focus the analysis on. Useful when you want to check saturation on a specific branch. |
| `-category <rate_category>` | Rate categories of interest. Relevant for models with gamma-distributed rate variations or FreeRate model. If the `-model` option includes rate variation (e.g., `+G4`), the `-category` should be a number between 1 and 4. |
| `-category_assignment` | Write assignment of the individual sites to the rate heterogeneity categories. |
| `-ufboot <number_of_replicates>` | Number of replicates for the ultrafast bootstrap analysis. Typically, a higher number like `1000` or `5000` is used. Ultrafast bootstrap provides rapid approximations to traditional bootstrap values. |
| `-boot <number_of_replicates>` | Number of replicates for traditional bootstrap analysis. This also computes a Maximum Likelihood (ML) tree and a consensus tree. Common value is `100`. |
| `-asr` | Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to a .asr.csv file. |
| `-output_suffix <output_suffix>` | Specify a suffix for the output file. |
| `-verbose` |  Enable verbose logging. |
| `-quiet`   |  Even no warnings on the the terminal. |
| `-dev`     |  Show the stack trace. |
