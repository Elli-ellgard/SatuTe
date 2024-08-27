# Examples Nucleotide Alignments

This directory provides examples to help users understand and utilize the SatuTe software tool effectively. The minimal input of SatuTe is a multiple sequence alignment, a model of sequence evolution with its parameters, and a phylogenetic tree. SatuTe parses these essential pieces of information from the output of [IQ-Tree](http://www.iqtree.org/). Therefore, the easiest way to get started with SatuTe is to use the option `-dir /path/to/iqtree_output`.  We provide here two examples subfolders: `dir_ML_tree` and `dir_true_tree`.

## Dataset Descriptions

### 1. Data

The `data` folder contains the raw data used for the constructed directory examples to demonstrate the functionality of SatuTe. Using a balanced 16-taxon trees, provided in `true_simulation_tree.tree`, on the taxon set A1, ..., A8, B1, ..., B8, we simulated a nucleotide alignment of length 1000bp under the JC model using [Seq-Gen](http://tree.bio.ed.ac.uk/software/seqgen/). The simulated alignment is stored in `example_dna.fasta`.  Furthermore the fully labelled tree is provided in `labelled_true_tree.tree`.

### 2. True Tree

The `dir_true_tree` folder includes  the output files of an [IQ-Tree](http://www.iqtree.org/) run using the following command:

```bash
iqtree2 -s example_dna.fasta -te labelled_true_tree.tree -m JC -blfix
```

This scenario represents the statistically correct situation: the phylogenetic signal of an simulated alignment is tested. The alignment was simulated along the true (simulation) tree with its branch lengths.

### 3. Maximum Likelihood Tree

The `dir_ML_tree` folder includes  the output files of an [IQ-Tree](http://www.iqtree.org/) run using the following command:

```bash
iqtree2 -s example_dna.fasta
```

With this command, the  maximum likelihood (ML) tree from a sequence alignment (example_dna.fasta) with the best-fit model automatically selected by ModelFinder will be inferred, a typical procedure in phylogenomics.

## Getting Started with SatuTe

To use the software tool with the provided iqtree runs, follow these steps:

1. **True Tree**:
    - To evaluates whether the alignment provides enough phylogenetic information to justify the branches in the tree, run SatuTe with the following command:

        ```bash
        satute -dir ./dir_true_tree
        ```

    - The SatuTe output files are directly written into the directory. There will be a log file `satute.log`, a test results file `satute.csv` containing the results for all branches, and a nexus file `satute.nex` containing the results as metadata. Additionally, there is a file `satute.components.csv` that outputs the components for each test statistic. For more information, refer to the SatuTe manual. For comparison the files are provided in the subfolder `SatuTe_results`.
    - We consider now the output file `satute.csv`. This file provides a comprehensive overview of the saturation test results for a specific branch or all branches

        We get the following results:
        | branch                       | ... | **z_score** | p_value | ...|**decision_test** | ... | branch_length | number_of_sites |
        |----------------------------|-----|-------------|---------|---|-------------------|-----|---------------|-----------------|
        | ...                        | ... | ...         | ...   | ...  | ...               | ... | ...           | ...             |
        | (NodeA5678*, NodeAroot*)   | ... | **38.7316** | 0   | ...    | **Informative**   | ... | 0.05108455    | 1000            |
        | (NodeBroot*, NodeAroot*)   | ... | **1.6294**  | 0.0516  | ... | **Saturated**     | ... | 2             | 1000            |
        | (NodeA56*, NodeA5678*)     | ... | **45.1553** | 0       | ... |**Informative**   | ... | 0.07450149    | 1000            |
        | ...                        | ... | ...         | ...   | ...  | ...               | ... | ...           | ...             |

        Only the central  branch AB is saturated.

2. **ML Tree**:

     - By default, the significance level is set to 5%. Analogously to the first example, run SatuTe adjusting the significance level to a more stringent one:

         ```bash
         satute -dir ./dir_ML_tree -alpha 0.01
         ```

     - The SatuTe output files are again directly written into the directory.  For comparison the files are provided in the subfolder `SatuTe_results`.
     - Reconstructing a ML tree from an alignment and then determining which branches of this ML tree are supported by the same alignment is a circular analysis that leads to an inflation of type I error. Therefore, a Bonferroni correction is necessary. Now, we consider the following columns of the output file `satute.csv`.

        | branch                       | ... | **z_score** | ... | **z_alpha_bonferroni_corrected** | **decision_bonferroni_corrected** | ... | branch_length | number_of_sites |
        |----------------------------|-----|-------------|-----|-----------------------|----------------------------------|-----|---------------|-----------------|
        | ...                        | ... | ...         | ... | ...                   | ...                              | ... | ...           | ...             |
        | (Node5*, Node4*)           | ... | **38.7977** | ... | **3.5293**            | **Informative**                  | ... | 0.0043852088  | 1000            |
        | (Node8*, Node4*)           | ... | **2.5172**  | ... | **3.6047**            | **Saturated**                    | ... | 2.2904773271  | 1000            |
        | (Node6*, Node5*)           | ... | **45.2164** | ... | **3.384**             | **Informative**                  | ... | 0.0702957449  | 1000            |
        | ...                        | ... | ...         | ... | ...                   | ...                      | ... | ...           | ...             |
  
        Again, the central branch is saturated.

## Further Commands

Given a path to an IQ-Tree executable, SatuTe runs IQ-Tree with default options to generate the required data for the test of saturation, namely  a multiple sequence alignment (MSA), a model of sequence evolution with its parameters, and a phylogenetic tree. The following commands can be used to generate the same results using other SatuTe options.

### True Tree

In this case, use the SatuTe mode  MSA+Model+Tree  with additional IQ-Tree arguments:

```bash
satute -msa ./data/example_dna.fasta -model JC -tree labelled_true_tree.tree \
        -iqtree iqtree2 -add_iqtree_options "-blfix"
```

### ML Tree

Here, use the SatuTe mode MSA:

```bash
satute -msa ./data/example_dna.fasta -iqtree path_to_iqtree_exe -alpha 0.01
```
