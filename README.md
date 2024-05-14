# Satute: Branch Saturation Testing Tool

## Introduction
Introduction
Satute is a specialized tool developed to measure the phylogenetic information shared by subtrees within a phylogeny, enabling the detection of branch saturation. In phylogenetic reconstruction, assessing the reliability of inferred trees and the data they are derived from is crucial. Saturation occurs when multiple substitutions obscure true genetic distances, potentially leading to artifacts and errors in phylogenetic analyses.

Satute addresses this by implementing a formal measure of phylogenetic information, facilitating a test for branch saturation. It wraps around the IQ-TREE software, leveraging its capabilities for likelihood calculations and evolutionary modeling. Users can input a Multiple Sequence Alignment (MSA), an evolutionary model, and a phylogenetic tree. Satute then performs its analysis, providing insights into the reliability of specific branches within the tree, thus helping researchers evaluate the robustness of their phylogenetic reconstructions.

By employing Satute, researchers can detect and quantify saturation, making informed decisions about the accuracy and stability of their phylogenetic trees.
## Usage: **Basic Usage**

1. **Using a Directory**:
   If you've previously run IQ-TREE and have a directory with the output files, you can provide the directory using the `-dir` option. This way, Satute will use the existing output without needing to rerun IQ-TREE: `bash python satute_cli.py -dir /path/to/iqtree/output/`
2. **Using a Multiple Sequence Alignment (MSA)**:
   As an alternative, you can provide a multiple sequence alignment (`-msa`) and the path to IQ-TREE (`-iqtree`). Then the best-fit evolutionary model will be identified using Modelfinder (as inmplemented in IQ-Tree) and a maximum likelihood tree will be inferred. IQ-Tree will run only with necessary options. For specific option choices, please run IQ-Tree separately and use the option `-dir` afterwards. Furthermore you are able to secify the tree (`-tree`) and the model of evolution (`-model`) togetheer with the option `-msa`, a typical command could look like: `bash python satute_cli.py -msa ./test/cassius/toyExample.phy -tree ./test/cassius/toyExample.phy.treefile -model GTR+G4 -iqtree iqtree`
   **Advanced Usage**
3. **Bootstrap Analysis**

   Ultrafast bootstrap analysis can be run using the `-ufboot` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -ufboot 1000
   ```

   Traditional bootstrap analysis can be performed using the `-boot` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -boot 100
   ```

   2.**Specifying an Edge for Analysis**
   If you want to focus the analysis on a specific branch or edge, use the `-edge` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -edge "(Node1, Node2)"
   ```

## Potential Errors and Warnings

1. **InvalidDirectoryError**: Thrown when the provided directory either does not exist or is empty. Ensure the directory path is correct and contains necessary IQ-TREE output files.

2. **NoAlignmentFileError**: Indicates that no multiple sequence alignment file was found in the specified directory. Ensure your directory contains the MSA file.

3. **ValueError**: Can occur in multiple scenarios:
   - If only the `-msa` and `-tree` options are used without specifying a model.

## Invalid Command Combinations

Certain combinations of command-line arguments are invalid:

1. **Directory with Model, MSA, Tree, Ufboot, Boot**: Providing an input directory with `-dir` shouldn't be combined with specifying a msa, a model, a tree, ufboot or boot option.
2. **Model and Tree without MSA**: Just providing the `-model` and `-tree` without a msa (`-msa`) is insufficient.

3. **MSA+Model+Tree with ufboot or boot option**: In the msa+model+tree mode, the inference is not re-done again, such that no ufboot and boot values can be determined.

4. **Edge without MSA**: The `-edge` option, used to focus the analysis on a specific branch, requires at least the `-msa` option or `-dir` option.

---

### Usage of `-add_iqtree_options`

The `-add_iqtree_options` flag allows you to specify additional options for the IQ-TREE run if necessary. This provides flexibility to customize the IQ-TREE execution by including specific commands that are not covered by the predefined arguments. You can use multiple additional options with this flag, particularly in scenarios where you are using the MSA option alone, the MSA option combined with a model, or the MSA option combined with both a tree and a model.

Here are some examples of how you can use this flag:

- **Using MSA alone**:
  ```sh
  python satute_cli.py -msa /path/to/alignment.fasta -add_iqtree_options "-alninfo"
  ```
  This command will print alignment site statistics to a `.alninfo` file.

- **Using MSA with a model**:
  ```sh
  python satute_cli.py -msa /path/to/alignment.fasta -model GTR+G4 -add_iqtree_options "-blfix"
  ```
  This command will fix branch lengths of the tree passed via `-tree` or `-te`.

- **Using MSA with a tree and a model**:
  ```sh
  python satute_cli.py -msa /path/to/alignment.fasta -tree /path/to/treefile.tree -model HKY -add_iqtree_options "-blmin 0.00001 -blmax 5"
  ```
  In this command, several options are used:
  - `-blmin`: Specifies the minimum branch length (default is the smaller of 0.000001 and 0.1/alignment_length).
  - `-blmax`: Specifies the maximum branch length (default is 10).

Additionally, when a user provides a given tree file using the `-tree` option, IQ-TREE automatically performs estimations based on that tree, streamlining the process. This means that, with a provided tree, IQ-TREE will automatically estimate the likelihood and perform other necessary computations without requiring additional input from the user for these steps. If you do not want IQ-TREE to estimate branch lengths automatically, you can use the `-blfix` option within `-add_iqtree_options` to fix the branch lengths.

Here are some useful additional IQ-TREE options that can be used with the `-add_iqtree_options` flag:
- `-alninfo`: Print alignment site statistics to a `.alninfo` file.
- `-blfix`: Fix branch lengths of the tree passed via `-tree` or `-te`.
- `-blmin`: Specify the minimum branch length (default: the smaller of 0.000001 and 0.1/alignment_length).
- `-blmax`: Specify the maximum branch length (default: 10).


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

---
