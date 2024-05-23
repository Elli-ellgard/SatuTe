# Satute: Branch Saturation Testing Tool

## Introduction
<div style="text-align: center;">
  <img src="./docs//figure_1_2024_05_22-1.png" alt="Example Image" />
  <p><em><p>
  <em>Figure 1: In the alignment of sequences from different species, each column represents a pattern. For a given tree, branch AB splits the tree into subtrees TA and TB, dividing each pattern into subpatterns. These subpatterns are used to compute likelihood vectors L(&#8706;A) and L(&#8706;B), and the scalar product <span style="font-family: 'Times New Roman', Times, serif;">C<sub>1</sub><sup>&#8706;</sup></span>. The average scalar product <span style="font-family: 'Times New Roman', Times, serif;">&#770;C<sub>1</sub></span> is computed, and using its variance and the number of sites, the z-score Z is calculated. Satute uses the z-score to test if branches are informative or saturated, based on a chosen significance level.</em>
</p>
</em><p></div>

Satute is a specialized tool designed to evaluate the phylogenetic information shared by subtrees within a phylogeny, enabling the detection of branch saturation. Assessing the reliability of inferred phylogenetic trees and their underlying data is crucial in phylogenetic reconstruction. Satute facilitates this by implementing a formal measure of phylogenetic information to test for branch saturation. Currently, users can either input a finished IQ-TREE run or use IQ-TREE on a Multiple Sequence Alignment (MSA) to apply the saturation test. Satute performs its analysis, providing insights into the reliability of specific branches within the tree, helping researchers assess the robustness of their phylogenetic reconstructions. By employing Satute, researchers can detect and quantify saturation, making informed decisions about the accuracy and stability of their phylogenetic trees. For more details on IQ-TREE, visit [IQ-TREE](http://www.iqtree.org) [Nguyen et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4271533/).

## Getting started: **Installing Satute**

To install Satute, you can use pipx, the Python package installer. Follow the instructions below to install the tool on your system.

pipx is a tool that allows you to install and run Python applications in isolated environments, ensuring each application and its dependencies are kept separate to avoid conflicts.

### Steps to Install Satute using pipx:

1. **Install Satute using pipx:**
   Once pipx is installed, you can use it to install Satute:
   ```sh
   pipx install satute
   ```
For more detailed instructions and information about pipx, refer to the [official pipx documentation](https://pipxproject.github.io/pipx/).

Using pipx ensures that Satute and its dependencies are installed in an isolated environment, minimizing potential conflicts with other Python packages on your system.
### Prerequisites

* Python 3.6 or higher
* `pipx` (Python package installer)

### Verifying the Installation

After the installation is complete, you can verify that Satute has been installed correctly by running the following command:

```bash
  satute version
```

You should see the version number of Satute printed on the screen, confirming that the installation was successful or to see if the installation has been successful you can also try to call.

```bash
  satute -h
```

There you will see all the options that Satute is providing.

## Usage: **Basic Usage**

Satute provides several ways to apply the saturation test on a given multiple sequence alignment, tree, and model. The easiest and most recommended approach is to use it directly on a finished IQ-TREE run. For details please go to the official website of IQ-Tree: http://www.iqtree.org/

1. **Using a Directory**:
If you've previously run IQ-TREE and have a directory with the output files, you can provide the directory using the `-dir` option. This way, Satute will use the existing output without needing to rerun IQ-TREE. However, if the tree has been inferred with a model with rate categories and the `-wspr` option has not been used, further analyses cannot be performed because the `.siteprob` file will not be created, making the assignment of each site to rate categories impossible. In such cases, we suggest rerunning IQ-TREE with the `-wspr` option.

   ```bash
   satute -dir ./examples/example_data_dna_iqtree_run/
   ```

When the run is completed, you will find several files in the output directory. Satute generates a `satute.log` file that logs the details of the run. If no specific rate category is defined, the following files will be created for the entire alignment and the given model:

* `satute.csv`
* `satute.components.csv`
* `satute.nex`

If a specific rate category is defined, these three files will be created for each rate category:

* `satute.<rate_category>.csv`
* `satute.<rate_category>.components.csv`
* `satute.<rate_category>.nex`

## satute.log
The `satute.log` file provides a comprehensive record of the steps and processes performed by the Satute tool during its execution. It includes details on the initialization and configuration, the substitution model used, spectral decomposition results, and the analysis execution. Additionally, it logs the writing of results to various output files and provides a summary of the number of sites corresponding to each rate category, ensuring a transparent and traceable analysis process.
## satute.csv file
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

## satute.components file

### Components File Description

| Column Name     | Description                                                       |
| --------------- | ----------------------------------------------------------------- |
| Edge            | The branch or edge in the tree being analyzed                     |
| coefficient     | The coefficient value for the site in the specified rate category |
| sample_variance | The variance of the coefficient for the site                      |
| rate            | The rate category                                                 |
| site            | The specific site in the alignment being analyzed                 |

### Description

The components file contains the variance and the coherence values for each site in the alignment for a specific edge in the tree. Each row represents a site with its corresponding coefficient, variance, and rate category for the edge "(t7, Node1*)".

## .satute.nex file

### Description of the NEXUS File

The NEXUS file contains two main sections: `TAXA` and `TREES`.

#### TAXA Section

Lists the 7 taxa included in the analysis:

```
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
```

### TREES Section

Defines the phylogenetic tree with metadata for each branch:

```
BEGIN TREES;
Tree tree1 = (t7:2.16965e-06[&coefficient_value=0.0889,standard_error_of_mean=0.1493,test_statistic=0.595,p_value=0.2759,z_alpha=1.6449,decision_test=Saturated,z_alpha_corrected=2.394,decision_corrected_test_tips=Saturated,decision_test_tip2tip=SatuT2T,branch_length=2.1696521999999996e-06,number_of_sites=8,rate_category=p4],...);
END;
```

Each branch includes:

* Coefficient value, standard error, test statistic, and p-value.
* Decisions based on statistical tests (e.g., Informative, Saturated).
* Branch length, number of sites, and rate category.

## Other Input Types to Use Satute

1. **Using a Multiple Sequence Alignment (MSA)**:
As an alternative, you can provide a multiple sequence alignment (`-msa`) and the path to IQ-TREE (`-iqtree`). If no model is specified, ModelFinder, which is included in IQ-TREE [Kalyaanamoorthy et al., 2017], will infer the best-fitting evolutionary model to conduct the analysis. IQ-TREE will run with only the necessary options. For specific option choices, please run IQ-TREE separately and use the `-dir` option afterwards. For the inference of the tree and other computations, `-iqtree` will be needed. By default, the path to IQ-TREE is set to `iqtree2`. If IQ-TREE is located elsewhere, you can specify its path using the `-iqtree` flag.

   **Example without specifying a model**:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -iqtree iqtree
   ```

   **Example specifying a model**:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -iqtree iqtree
   ```

1. **Using a MSA+Model+Tree**:
You can also specify a tree file (`-tree`) along with the multiple sequence alignment (`-msa`). However, without specifying a model (`-model`), this will lead to an error. Ensure that the model is specified when providing a tree file. Additionally, note that branch lengths will be reestimated during the analysis. If you want to fix the branch lengths, you need to add the `-blfix` option using the `-add_iqtree_options` flag. Below are example commands:

   **Example with branch length reestimation**:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -tree ./test/cassius/toyExample.phy.treefile -model GTR+G4 -iqtree iqtree
   ```

   **Example fixing branch lengths**:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -tree ./test/cassius/toyExample.phy.treefile -model GTR+G4 -iqtree iqtree -add_iqtree_options "-blfix"
   ```

### Usage of `-add_iqtree_options`

The `-add_iqtree_options` flag allows you to specify additional options for the IQ-TREE run if necessary. This provides flexibility to customize the IQ-TREE execution by including specific commands that are not covered by the predefined arguments. You can use multiple additional options with this flag, particularly in scenarios where you are using the MSA option alone, the MSA option combined with a model, or the MSA option combined with both a tree and a model.

Here are some examples of how you can use this flag:

**Using MSA alone**:

  ```sh
  satute -msa /path/to/alignment.fasta -add_iqtree_options "-alninfo"
  ```
  This command will print alignment site statistics to a `.alninfo` file.

**Using MSA with a model**:
  
  ```sh
  satute.py -msa /path/to/alignment.fasta -model GTR+G4 -add_iqtree_options "-blfix"
  ```

  This command will fix branch lengths of the tree passed via `-tree` or `-te`.

**Using MSA with a tree and a model**:
  
  ```sh 
  satute -msa /path/to/alignment.fasta -tree /path/to/treefile.tree -model HKY -add_iqtree_options "-blmin 0.00001 -blmax 5"
  ```

  In this command, several options are used:
  - `-blmin`: Specifies the minimum branch length (default is the smaller of 0.000001 and 0.1/alignment_length).
  - `-blmax`: Specifies the maximum branch length (default is 10).

* **Advanced Usage**

1. **Bootstrap Analysis**

   Ultrafast bootstrap analysis can be run using the `-ufboot` option:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -ufboot 1000
   ```

   Traditional bootstrap analysis can be performed using the `-boot` option:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -boot 100
   ```

   2.**Specifying an Edge for Analysis**
   If you want to focus the analysis on a specific branch or edge, use the `-edge` option:

   ```bash
   satute -msa ./test/cassius/toyExample.phy -model GTR+G4 -edge "(Node1, Node2)"
   ```

## Potential Errors and Warnings

1. **InvalidDirectoryError**: Thrown when the provided directory either does not exist or is empty. Ensure the directory path is correct and contains necessary IQ-TREE output files.

2. **NoAlignmentFileError**: Indicates that no multiple sequence alignment file was found in the specified directory. Ensure your directory contains the MSA file.

3. **ValueError**: Can occur in multiple scenarios:
   * If only the `-msa` and `-tree` options are used without specifying a model.

## Invalid Command Combinations

Certain combinations of command-line arguments are invalid:

1. **Directory with Model, MSA, Tree, Ufboot, Boot**: Providing an input directory with `-dir` shouldn't be combined with specifying a msa, a model, a tree, ufboot or boot option.
2. **Model and Tree without MSA**: Just providing the `-model` and `-tree` without a msa (`-msa`) is insufficient.

3. **MSA+Model+Tree with ufboot or boot option**: In the msa+model+tree mode, the inference is not re-done again, such that no ufboot and boot values can be determined.

4. **Edge without MSA**: The `-edge` option, used to focus the analysis on a specific branch, requires at least the `-msa` option or `-dir` option.

### Command Line Arguments

This script accepts the following command-line arguments:

- `-dir <directory_path>`:
  - **Description**: Path to the input directory containing IQ-TREE output files. Use this option when you've already run IQ-TREE and want to avoid rerunning it. The directory should contain essential IQ-TREE output files including the .iqtree file, tree file(s), and possibly a .siteprob file.
  - **Type**: Valid directory
  - **Example**: `-dir /path/to/iqtree/output`

- `-tree <tree_file_path>`:

  - **Description**: Path to the input tree file in Newick or Nexus format. This tree will be used as the basis for the saturation analysis.
  - **Type**: Valid file
  - **Example**: `-tree /path/to/treefile.tree`

- `-msa <msa_file_path>`:

  - **Description**: Path to the Multiple Sequence Alignment (MSA) file you wish to analyze. The MSA can be in FASTA, NEXUS, PHYLIP, or TXT format.
  - **Type**: Valid file
  - **Example**: `-msa /path/to/alignment.fasta`

- `-iqtree <iqtree_path>`:

  - **Description**: Specifies the path to the IQ-TREE executable. If IQ-TREE is installed system-wide, just providing the executable name (`iqtree` or `iqtree2`) will suffice. Otherwise, give the complete path.
  - **Default**: `iqtree2`
  - **Type**: Path
  - **Example**: `-iqtree /usr/local/bin/iqtree2`

- `-model <evolution_model>`:

  - **Description**: Indicates the model of sequence evolution. Common models include `GTR`, `HKY`, etc. You can also specify rate heterogeneity and other model extensions, like `+G4` for gamma-distributed rates.
  - **Type**: String
  - **Example**: `-model GTR+G4`

- `-category <rate_category>`:

  - **Description**: Rate categories of interest. Relevant for models with gamma-distributed rate variations or FreeRate model. If the `-model` option includes rate variation (e.g., `+G4`), the `-category` should be a number between 1 and 4.
  - **Type**: Integer
  - **Example**: `-category 4`

- `-ufboot <number_of_replicates>`:

  - **Description**: Number of replicates for the ultrafast bootstrap analysis. Typically, a higher number like `1000` or `5000` is used. Ultrafast bootstrap provides rapid approximations to traditional bootstrap values.
  - **Type**: Integer
  - **Example**: `-ufboot 1000`

- `-boot <number_of_replicates>`:

  - **Description**: Number of replicates for traditional bootstrap analysis. This also computes a Maximum Likelihood (ML) tree and a consensus tree. Common values are `1000` or `5000`.
  - **Type**: Integer
  - **Example**: `-boot 1000`

- `-alpha <significance_level>`:

  - **Description**: Significance level for the saturation test. A common threshold is `0.05`, indicating a 5% significance level. Lower values make the test more stringent.
  - **Type**: Float
  - **Default**: `0.05`
  - **Example**: `-alpha 0.01`

- `-edge <edge_name>`:

  - **Description**: Specify a branch or edge name to focus the analysis on. Useful when you want to check saturation on a specific branch.
  - **Type**: String
  - **Example**: `-edge branch1`

- `-output_suffix <output_suffix>`:

  - **Description**: Specify a suffix for the output file.
  - **Type**: String
  - **Default**: `""`
  - **Example**: `-output_suffix _analysis`

- `-add_iqtree_options <additional_option>`:

  - **Description**: Specify additional options for the IQ-Tree run, if necessary.
  - **Type**: String
  - **Example**: `-add_iqtree_options "-nt AUTO"`

- `-asr`:

  - **Description**: Write ancestral sequences (by empirical Bayesian method) for all nodes of the tree to a .asr.csv file.
  - **Type**: Flag
  - **Example**: `-asr`

- `-category_assignment`:

  - **Description**: Write assignment of the individual sites to the rate heterogeneity categories.
  - **Type**: Flag
  - **Example**: `-category_assignment`

- `-verbose`:
  - **Description**: Enable verbose logging.
  - **Type**: Flag
  - **Example**: `-verbose`

