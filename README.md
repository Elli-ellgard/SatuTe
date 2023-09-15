
---

# Satute: Branch Saturation Testing Tool

## Introduction

Satute is a specialized tool designed to assess the saturation of branches in phylogenetic trees. By analyzing the substitution patterns across branches, Satute provides insights into potential artifacts and errors in the tree, helping researchers make more informed decisions about the reliability of their phylogenetic reconstructions.

Satute wraps around the IQ-TREE tool and adds functionalities specific to branch saturation testing. The user can provide a Multiple Sequence Alignment (MSA), a tree, and specify an evolutionary model. Satute will then perform its analysis using this input, potentially leveraging IQ-TREE for some computations.

## Usage:

**Basic Usage**

To run Satute with the most common settings, you can provide the multiple sequence alignment (`-msa`), the tree (`-tree`), the model of evolution (`-model`), and the path to IQ-TREE (`-iqtree`):

```bash
python satute_cli.py -msa ./test/cassius/toyExample.phy -tree ./test/cassius/toyExample.phy.treefile -model GTR+G4 -iqtree iqtree
```

**Advanced Usage**

1. **Using a Directory**
   If you've previously run IQ-TREE and have a directory with the output files, you can provide the directory using the `-dir` option. This way, Satute will use the existing output without needing to rerun IQ-TREE:

   ```bash
   python satute_cli.py -dir /path/to/iqtree/output/
   ```


2. **Bootstrap Analysis**

   Ultrafast bootstrap analysis can be run using the `-ufboot` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -ufboot 1000
   ```

   Traditional bootstrap analysis can be performed using the `-boot` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -boot 1000
   ```

3. **Specifying an Edge for Analysis**

   If you want to focus the analysis on a specific branch or edge, use the `-edge` option:

   ```bash
   python satute_cli.py -msa ./test/cassius/toyExample.phy -model GTR+G4 -edge "(Node1, Node2)"
   ```

## Potential Errors and Warnings:

1. **InvalidDirectoryError**: Thrown when the provided directory either does not exist or is empty. Ensure the directory path is correct and contains necessary IQ-TREE output files.
  
2. **NoAlignmentFileError**: Indicates that no multiple sequence alignment file was found in the specified directory. Ensure your directory contains the MSA file.

3. **ValueError**: Can occur in multiple scenarios:
    - If the number of rate categories provided (`-nr`) doesn't match the number in the model.
    - If only the `-msa` and `-tree` options are used without specifying a model.

## Invalid Command Combinations:

Certain combinations of command-line arguments are invalid:

1. **Directory with Model**: Providing an input directory with `-dir` shouldn't be combined with specifying a model using `-model`.
   
2. **MSA without Model or Directory**: Just providing the `-msa` without either a model (`-model`) or directory (`-dir`) is insufficient.

3. **Edge without MSA**: The `-edge` option, used to focus the analysis on a specific branch, requires the `-msa` option.

---
### Arguments:

- **-dir**:  
  **Type**: Path  
  **Default**: None  
  **Description**: Specifies the path to an existing directory containing IQ-TREE output files. This is useful when you've already run IQ-TREE and would like to avoid rerunning it. The directory should contain essential IQ-TREE output files including the `.iqtree` file, tree file(s), and, if the model uses multiple rate categories, a `.siteprob` file.  
  **Example**: `-dir /path/to/iqtree_output/`

- **-tree**:  
  **Type**: Path  
  **Default**: None  
  **Description**: Provides the path to an input tree file. This tree will be used for the saturation analysis. The tree file can be in Newick or Nexus format.  
  **Example**: `-tree /path/to/treefile.nwk`

- **-msa**:  
  **Type**: Path  
  **Default**: None  
  **Description**: Path to the Multiple Sequence Alignment (MSA) file that you want to analyze. This file can be in FASTA, NEXUS, PHYLIP, or TXT format.  
  **Example**: `-msa /path/to/alignment.fasta`

- **-iqtree**:  
  **Type**: Path  
  **Default**: None  
  **Description**: Specifies the path to the IQ-TREE executable. If you've installed IQ-TREE system-wide, just the executable name (`iqtree` or `iqtree2`) will suffice. Otherwise, provide the full path.  
  **Example**: `-iqtree /path/to/iqtree2`

- **-model**:  
  **Type**: String  
  **Default**: None  
  **Description**: Represents the model of sequence evolution to be used. Common models include `GTR`, `HKY`, and more. Rate heterogeneity and other model extensions can be appended, such as `+G4` for gamma-distributed rates.  
  **Example**: `-model GTR+G4`

- **-nr**:  
  **Type**: Integer  
  **Default**: None  
  **Description**: Specifies the number of rate categories for the model. This is especially relevant for models with gamma-distributed rate variations. If the `-model` option includes a rate variation (e.g., `+G4`), the `-nr` should match the number specified in the model.  
  **Example**: `-nr 4`

- **-ufboot**:  
  **Type**: Integer (>=1000)  
  **Default**: None  
  **Description**: Indicates the number of replicates to be used for ultrafast bootstrap analysis. Generally, a higher number of replicates provides more robust support values. Typically, values like `1000` or `5000` are used.  
  **Example**: `-ufboot 1000`

- **-boot**:  
  **Type**: Integer  
  **Default**: None  
  **Description**: Specifies the number of replicates for traditional bootstrap analysis. This option also triggers the computation of a Maximum Likelihood (ML) tree and a consensus tree. As with ultrafast bootstrap, values like `1000` or `5000` are commonly used.  
  **Example**: `-boot 1000`

- **-alpha**:  
  **Type**: Float  
  **Default**: 0.05  
  **Description**: Determines the significance level for the saturation test. A common value is `0.05`, representing a 5% significance level. Lower values make the test more stringent.  
  **Example**: `-alpha 0.05`

- **-edge**:  
  **Type**: String  
  **Default**: None  
  **Description**: If you want to focus the analysis on a specific branch or edge, you can specify its name using this option. This is particularly useful when you want to assess saturation on a particular branch of interest.  
  **Example**: `-edge EdgeName`

---
