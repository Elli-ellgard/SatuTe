# Files in the Folder

We have two folders with two datasets: one for Amino Acids (`example_satute_aa`) and one for DNA (`example_satute_dna`). Both datasets are derived from simulated data. The runs without the Gamma model have been copied, and an IQ-TREE run has been executed on them, allowing you to test Satute directly on a finished IQ-TREE run for both dataset types.

## Directory Structure

### Amino Acids Dataset (`example_satute_aa`)
- **Description:** Contains simulated alignment data for amino acids.
- **Subdirectories:**
  - **simtree:** Contains the original trees used for the alignment simulations with seqgen, both with and without labels. The trees are provided in Phylip format (`simtree.tree`) and labeled format (`simtree-labelled.tree`). These trees were used to generate the alignments with seq-gen under the LG model with the model's stationary frequencies.
  - **cbl3.5_A8+B8_1500bp_LG_iqtree_run:** Contains the IQ-TREE run output, ready for direct use with Satute. The README in this directory explains that the trees have two subtrees of 8 taxa each, named A1...A8 and B1...B8. The branches inside the subtrees are set to 0.1, and the internal branch linking the subtrees is set to 3.5 substitutions per site. The DNA alignments of length 1500bp were simulated using seq-gen under the LG model with the model's stationary frequencies.
  - **cbl3.5_A8+B8_1500bp_LG:** Contains the alignment data for testing with MSA, MSA+Model, and MSA+Model+Tree options.
  - **cbl3.5_A8+B8_1500bp_LG+G:** Contains the alignment data create with a Gamma Model


### DNA Dataset (`example_satute_dna`)
- **Description:** Contains simulated alignment data for DNA.
- **Subdirectories:**
  - **simtree:** Contains the original trees used for the alignment simulations with seqgen, both with and without labels. The trees are provided in Phylip format (`simtree.tree`) and labeled format (`simtree-labelled.tree`). These trees were used to generate the alignments with seq-gen under the JC model.
  - **tree_cbl3.5_A8+B8_1500bp_JC_iqtree_run:** Contains the IQ-TREE run output, ready for direct use with Satute. The README in this directory explains that the trees have two subtrees of 8 taxa each, named A1...A8 and B1...B8. The branches inside the subtrees are set to 0.1, and the internal branch linking the subtrees is set to 3.5 substitutions per site. The DNA alignments of length 1500bp were simulated using seq-gen under the JC model.
  - **tree_cbl3.5_A8+B8_1500bp_JC:** Contains the alignment data for testing with MSA, MSA+Model, and MSA+Model+Tree options.
  - **tree_cbl3.5_A8+B8_1500bp_JC+G:** Contains the alignment data to create with a Gamma Model.
  
## Usage Instructions

### Using Finished IQ-TREE Runs

To test Satute directly on a finished IQ-TREE run, use the following command:

```bash
satute -dir ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG_iqtree_run
```

### Using Multiple Sequence Alignment (MSA)

To use the datasets with MSA, MSA+Model, and MSA+Model+Tree options, use the following commands:

**Using MSA:**
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG
```

**Using MSA with a specified model:**
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG -model LG
```

**Using MSA with a specified model and tree:**
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG -model LG -tree ./examples_satute_aa/simtree/simtree.tree
```

## Trees Used for Simulation

In all directories, there is a `simtree` directory containing the trees used for the simulations. These trees are provided both with and without labels. The labeled trees (`simtree-labelled.tree`) include inner nodes labeled with "Node1*", which can be useful for specific analyses.

## Examples

### Amino Acids Dataset Example

Run Satute on the finished IQ-TREE run:
```bash
satute -dir ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG_iqtree_run
```

Use the MSA file directly:
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG
```

Specify the model:
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG -model LG
```

Specify the model and tree:
```bash
satute -s ./examples_satute_aa/cbl3.5_A8+B8_1500bp_LG -model LG -tree ./examples_satute_aa/simtree/simtree.tree
```

### DNA Dataset Example

Run Satute on the finished IQ-TREE run:
```bash
satute -dir ./examples_satute_dna/tree_cbl3.5_A8+B8_1500bp_JC_iqtree_run
```

Use the MSA file directly:
```bash
satute -s ./examples_satute_dna/tree_cbl3.5_A8+B8_1500bp_JC
```

Specify the model:
```bash
satute -s ./examples_satute_dna/tree_cbl3.5_A8+B8_1500bp_JC -model JC
```

Specify the model and tree:
```bash
satute -s ./examples_satute_dna/tree_cbl3.5_A8+B8_1500bp_JC -model JC -tree ./examples_satute_dna/simtree/simtree.tree
```

