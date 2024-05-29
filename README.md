# Software Manual for Satute

## Table of Contents

1. [Introduction](#1-introduction)
   - [Overview](#overview)
   - [Requirements](#requirements)
   - [Main Workflow](#main-workflow)
2. [Installation](#2-installation)
   - [Prerequisites](#prerequisites)
   - [Using Pipex](#install-satute-using-pipx)
   - [Verifying the Installation](#verifying-the-installation)
3. [Getting Started](#getting-started)
   - [Minimal command-line examples](#minimal-)
   - [Initial Setup](#initial-setup)
4. [Satute Output](#4-satute-output)
   - [Log file](#log-file)
   - [Test Results](#test-results)
   - [Components of Test Statistic](#components-of-test-statistic)
   - [Nexus file](#nexus-file)
5. [Additional Features](#5-additional-features)
   - [Bootstrap Analysis](#bootstrap-analysis)
   - [Edge Specification](#edge-specification)
   - [Category Features](#category-features)
   - [Ancestral sequence reconstruction](#ancestral-sequence-reconstruction)
   - [Others](#others)
6. [Troubleshooting](#6-troubleshooting)
   - [Potential Errors and Warnings](#potential-errors-and-warnings)
   - [Invalid Command Combinations](#invalid-command-combinations)
   - [Support and Contact Information](#support-and-contact-information)

## 1. Introduction

### Overview

Welcome to the Satute manual. This document provides comprehensive information on how to install, configure, and use Satute effectively.

Satute is a Python-based tool designed to test for branch saturation in phylogenetic analyses. Saturation occurs when multiple substitutions obscure true genetic distances, potentially leading to artifacts and errors. Assessing the reliability of inferred phylogenetic trees and the data they are derived from is crucial in phylogenetic reconstruction. The developed test for branch saturation measures the phylogenetic information shared by subtrees within a phylogeny, enabling the detection of branch saturation. By using Satute, researchers can detect and quantify saturation, thereby making informed decisions about the accuracy and stability of their phylogenetic trees.

### Requirements

The minimal input of Satute is a multiple sequence alignment, a model of sequence evolution with its parameters, and a phylogenetic tree. Satute parses these essential pieces of information from the output of [IQ-TREE](http://www.iqtree.org/), an efficient software for phylogenomic inference . While we strongly recommend running IQ-TREE separately with customised options, Satute can also use an IQ-TREE executable to generate any missing information using default settings.

**Technical Requirements:**

    1. Python: 3.6 or higher
    2. IQ-Tree: 2.2.2.3 or higher

### Main Workflow

The main function of Satute operates as follows: Given the required input, Satute first calculates the spectral decomposition of the rate matrix and determines the likelihood vectors for each node in the tree. It then performs the test for branch saturation on a user-selected branch or on every branch of the tree, as described in the relevant literature. The program outputs the test results and its components in different CSV files and a Nexus file (see section [Satute Output](#satute-output)).

In cases where a model of rate heterogeneity is used, Satute assigns each site to the rate category with the highest posterior probability. The alignment is then split by Satute. For each category, the test for branch saturation is employed on the rescaled phylogenetic tree given the subalignment.

## 2. Installation

Satute is available as a python package from pypi and can be normally installed via pip.
We recommend to use [pipx](https://pipx.pypa.io/stable/) to install Satute as a standalone command line tool. Using pipx ensures that Satute and its dependencies are installed in an isolated environment, minimizing potential conflicts with other Python packages on your system.

### Prerequisites

- Python 3.6 or higher
- `pipx` (Python package installer)

### Install Satute using pipx

1. **Install pipx:**  If you don't have pipx installed, you can install it using pip:

    ```sh
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

### Using a Directory:

   If you've previously run IQ-TREE and have a directory with the output files, you can provide the directory using the `-dir` option. This way, Satute will use the existing output without needing to rerun IQ-TREE: `bash python satute_cli.py -dir /path/to/iqtree/output/`
   
### **Using a Multiple Sequence Alignment (MSA)**:
   As an alternative, you can provide a multiple sequence alignment (`-msa`) and the path to IQ-TREE (`-iqtree`). Then the best-fit evolutionary model will be identified using Modelfinder (as inmplemented in IQ-Tree) and a maximum likelihood tree will be inferred. IQ-Tree will run only with necessary options. For specific option choices, please run IQ-Tree separately and use the option `-dir` afterwards. Furthermore you are able to secify the tree (`-tree`) and the model of evolution (`-model`) togetheer with the option `-msa`, a typical command could look like: `bash python satute_cli.py -msa ./test/cassius/toyExample.phy -tree ./test/cassius/toyExample.phy.treefile -model GTR+G4 -iqtree iqtree`

| Option | Description  | Example  |
| ----------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `-dir <directory_path>`                   | Path to the input directory containing IQ-TREE output files. Use this option when you've already run IQ-TREE and want to avoid rerunning it. The directory should contain essential IQ-TREE output files including the .iqtree file, tree file(s), and possibly a .siteprob file. | `-dir /path/to/iqtree/output`    |
| `-tree <tree_file_path>`                  | Path to the input tree file in Newick or Nexus format. This tree will be used as the basis for the saturation analysis.                                                                                                                                                           | `-tree /path/to/treefile.tree`   |
| `-msa <msa_file_path>`                    | Path to the Multiple Sequence Alignment (MSA) file you wish to analyze. The MSA can be in FASTA, NEXUS, PHYLIP, or TXT format.                                                                                                                                                    | `-msa /path/to/alignment.fasta`  |
| `-iqtree <iqtree_path>`                   | Specifies the path to the IQ-TREE executable. If IQ-TREE is installed system-wide, just providing the executable name (`iqtree` or `iqtree2`) will suffice. Otherwise, give the complete path.                                                                                    | `-iqtree /usr/local/bin/iqtree2` |
| `-model <evolution_model>`                | Indicates the model of sequence evolution. Common models include `GTR`, `HKY`, etc. You can also specify rate heterogeneity and other model extensions, like `+G4` for gamma-distributed rates.                                                                                   | `-model GTR+G4`                  |
| `-alpha <significance_level>`             | Significance level for the saturation test. A common threshold is `0.05`, indicating a 5% significance level. Lower values make the test more stringent.                                                                                                                          | `-alpha 0.01`                    |
| `-add_iqtree_options <additional_option>` | Specify additional options for the IQ-Tree run, if necessary.                                                                                                                                                                                                                     | `-add_iqtree_options "-nt AUTO"` |


## 4. Satute Output

TODO: table




## Log file

The `satute.log` file provides a comprehensive record of the steps and processes performed by the Satute tool during its execution. It includes details on the initialization and configuration, the substitution model used, spectral decomposition results, and the analysis execution. Additionally, it logs the writing of results to various output files and provides a summary of the number of sites corresponding to each rate category, ensuring a transparent and traceable analysis process.

## Test Results

In the satute.csv file one will find, for each edge the columns:

### Table Headers and Description

| Column Name                  | Description                                                    |
| ---------------------------- | -------------------------------------------------------------- |
| edge                         | The branch or edge in the tree being analyzed                  |
| coefficient_value            | The value of the coefficient calculated for the edge           |
| standard_error_of_mean       | The standard error of the mean for the coefficient             |
| test_statistic               | The test statistic value used to evaluate the edge             |
| p_value                      | The p-value indicating the significance of the test statistic  |
| z_alpha                      | The z-value corresponding to the alpha level for the test      |
| decision_test                | The decision based on the test statistic (e.g., Informative)   |
| z_alpha_corrected            | The corrected z-value considering multiple testing corrections |
| decision_corrected_test      | The decision based on the corrected z-value                    |
| decision_corrected_test_tips | The decision based on the corrected test for tips              |
| decision_test_tip2tip        | The decision based on the test for tip-to-tip comparisons      |
| branch_length                | The length of the branch or edge in the tree                   |
| number_of_sites              | The number of sites in the alignment associated with the edge  |
| rate_category                | The rate category for which the analysis was performed         |

## Components of Test Statistic

### Components File

| Column Name     | Description                                                       |
| --------------- | ----------------------------------------------------------------- |
| Edge            | The branch or edge in the tree being analyzed                     |
| coefficient     | The coefficient value for the site in the specified rate category |
| sample_variance | The variance of the coefficient for the site                      |
| rate            | The rate category                                                 |
| site            | The specific site in the alignment being analyzed                 |

### Description

The components file contains the variance and the coherence values for each site in the alignment for a specific edge in the tree. Each row represents a site with its corresponding coefficient, variance, and rate category for the edge "(t7, Node1*)".

## Nexus file

### Description of the NEXUS File

The NEXUS file contains two main sections: `TAXA` and `TREES`.

#### TAXA Section

Lists the 7 taxa included in the analysis:

```
#NEXUS
BEGIN TAXA;
    DIMENSIONS NTAX=7;
     TAXLABELS
        t7
        t3
        t2
        t5
        t6
        t1
        t4
    ;
END;

BEGIN TREES;
Tree tree1 = (t7:3.01328e-06,(((t3:1.55499,(t2:1.77629,t5:2.76104e-06)Node5*:0.377782[&p_value=0.0344,decision_test=Informative,decision_corrected_test_tips=Saturated])Node4*:0.368276,t6:2.16996e-06)Node3*:1.23617,t1:2.22639e-06)Node2*:1.05052,t4:1.85109)Node1*:0;
END;
```

Each branch includes:

* Decisions based on statistical tests p-value(e.g., Informative, Saturated).
* Branch length, number of sites, and rate category.


## 5. Additional Features

### Bootstrap Analysis

1. **Satute options:**

     | Option | Description  |
     | ------------ | ----------------------------------------------------------- |
     | `-ufboot <number_of_replicates>`          | Number of replicates for the ultrafast bootstrap analysis. Typically, a higher number like `1000` or `5000` is used. Ultrafast bootstrap provides rapid approximations to traditional bootstrap values. |
     | `-boot <number_of_replicates>`            | Number of replicates for traditional bootstrap analysis. This also computes a Maximum Likelihood (ML) tree and a consensus tree. Common values is `100`. |
     | |  |

2. **Examples:**

    Ultrafast bootstrap analysis can be run using the `-ufboot` option:

    ```bash
    satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -ufboot 1000
    ```

    Traditional bootstrap analysis can be performed using the `-boot` option:

    ```bash
    satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -boot 100
    ```

### Edge Specification

1. **Satute options:**

     | Option | Description  |
     | ------------ | ----------------------------------------------------------- |
     | `-edge <edge_name>`                       | Specify a branch or edge name to focus the analysis on. Useful when you want to check saturation on a specific branch.|
     | |  |

2. **Example:**
  
    If you want to focus the analysis on a specific branch or edge, use the `-edge` option:

    ```bash
    satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -edge "(Node1*, Node2*)"
    ```

### Category Features

1. **Satute options:**

     | Option | Description  |
     | ------------ | ----------------------------------------------------------- |
     | `-category <rate_category>`               | Rate categories of interest. Relevant for models with gamma-distributed rate variations or FreeRate model. If the `-model` option includes rate variation (e.g., `+G4`), the `-category` should be a number between 1 and 4.|
     | `-category_assignment`                    | Write assignment of the individual sites to the rate heterogeneity categories.|
     | |  |

2. **Example:**

     ....



### Ancestral sequence reconstruction

1. **Satute options:**

     | Option | Description  |
     | ------------ | ----------------------------------------------------------- |
     | `-asr`                                    | Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to a .asr.csv file. |
     | |  |

2. **Example:**


    The `-asr` option in Satute allows users to write ancestral sequences for all nodes of the tree to a `.asr.csv` file using the empirical Bayesian method. The `.asr.csv` file contains the posterior distributions of ancestral sequences for both the left and right nodes of the split trees at each edge. This feature is useful to gain more insight into the likelihoods of nodes that are separated by the edge being analyzed:

    ```bash
    satute -msa /path/to/alignment.fasta -model GTR+G4 -iqtree /usr/local/bin/iqtree2 -asr
    ```

    In this example, Satute will perform the analysis on the given multiple sequence alignment (`alignment.fasta`) using the specified evolutionary model (`GTR+G4`) and the IQ-TREE executable (`/usr/local/bin/iqtree2`). The ancestral sequences will be inferred and saved to a `.asr.csv` file.

3. **Result file `.asr.csv`:**

    The `.asr.csv` file contains the posterior distributions of ancestral sequences for both the left and right nodes of the split trees at each edge. The columns in the file are:

    | Column Name  | Description                                                             |
    | ------------ | ----------------------------------------------------------------------- |
    | `pA_left`    | Probability of nucleotide A in the left subtree at the specified site.  |
    | `pC_left`    | Probability of nucleotide C in the left subtree at the specified site.  |
    | `pG_left`    | Probability of nucleotide G in the left subtree at the specified site.  |
    | `pT_left`    | Probability of nucleotide T in the left subtree at the specified site.  |
    | `Node_left`  | Name of the node in the left subtree.                                   |
    | `pA_right`   | Probability of nucleotide A in the right subtree at the specified site. |
    | `pC_right`   | Probability of nucleotide C in the right subtree at the specified site. |
    | `pG_right`   | Probability of nucleotide G in the right subtree at the specified site. |
    | `pT_right`   | Probability of nucleotide T in the right subtree at the specified site. |
    | `Node_right` | Name of the node in the right subtree.                                  |
    | `Site`       | Site number in the alignment.                                           |
    | `Edge`       | Edge in the phylogenetic tree that splits the left and right subtrees.  |
    |||

    These columns represent the posterior distributions for the nodes on both sides of each edge in the phylogenetic tree. The probabilities are calculated for each site in the alignment, providing a detailed view of the ancestral sequence distributions at every edge split.

### Others

1. **Satute options:**

     | Option | Description  |
     | ------------ | ----------------------------------------------------------- |
     | `-output_suffix <output_suffix>`          | Specify a suffix for the output file.  |
     | `-verbose`                                | Enable verbose logging.| 
     | |  |

2. **Example:**

     ....


## 6. Troubleshooting

### Potential Errors and Warnings

1. **InvalidDirectoryError:**

    Thrown when the provided directory either does not exist or is empty. Ensure the directory path is correct and contains necessary IQ-TREE output files.

2. **NoAlignmentFileError:**

    Indicates that no multiple sequence alignment file was found in the specified directory. Ensure your directory contains the MSA file.

3. **ValueError:**

    Can occur in multiple scenarios:
    - If only the `-msa` and `-tree` options are used without specifying a model.

### Invalid Command Combinations

Certain combinations of command-line arguments are invalid:

1. **Directory with Model, MSA, Tree, Ufboot, Boot**: Providing an input directory with `-dir` shouldn't be combined with specifying a msa, a model, a tree, ufboot or boot option.

2. **Model and Tree without MSA**: Just providing the `-model` and `-tree` without a msa (`-msa`) is insufficient.

3. **MSA+Model+Tree with ufboot or boot option**: In the msa+model+tree mode, the inference is not re-done again, such that no ufboot and boot values can be determined.

4. **Edge without MSA or DIR**: The `-edge` option, used to focus the analysis on a specific branch, requires the `-msa` option or `-dir` option.


### Support and Contact Information

......