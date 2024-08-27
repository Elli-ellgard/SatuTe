# Example Amino Acids Alignments

This directory provides examples to help users understand and utilize the SatuTe software tool effectively. The minimal input of SatuTe is a multiple sequence alignment, a model of sequence evolution with its parameters, and a phylogenetic tree. SatuTe parses these essential pieces of information from the output of [IQ-Tree](http://www.iqtree.org/). Therefore, the easiest way to get started with SatuTe is to use the option `-dir /path/to/iqtree_output`. We provide here the subfolder: `dir_ML_tree`.

## Dataset Descriptions

### 1. Data

The `data` folder contains the raw data used for the constructed directory examples to demonstrate the functionality of SatuTe. Using a balanced 16-taxon trees, provided in `true_tree.tree`, on the taxon set A1, ..., A8, B1, ..., B8, we simulated amino acid alignment of length 1000bp under the LG+G model using [Seq-Gen](http://tree.bio.ed.ac.uk/software/seqgen/). The simulated alignment is stored in `example_aa.fasta`.  Furthermore the fully labelled tree is provided in `labelled_true_tree.tree`.

### 2. Maximum Likelihood Tree

The `dir_ML_tree` folder includes the output files of an [IQ-Tree](http://www.iqtree.org/) run using the following command:

```bash
iqtree2 -s example_aa.fasta -m LG+G -wspr
```

This command infers the maximum likelihood (ML) tree from a sequence alignment (example_aa.fasta), utilizing the best-fit model as automatically selected by ModelFinder, a standard approach in phylogenomic analysis.

## Getting Started with SatuTe

To use the software tool with the provided iqtree run, follow these steps:

1. **ML Tree**:

     - By default, the significance level is set to 5%. Analogously to the first example, run SatuTe adjusting the significance level to a more stringent one:

         ```bash
         satute -dir ./dir_ML_tree -alpha 0.01
         ```

     - The SatuTe output files for each rate category are directly written into the directory. For comparison the files are provided in the subfolder `SatuTe_results`.
     - Reconstructing a ML tree from an alignment and then determining which branches of this ML tree are supported by the same alignment is a circular analysis that leads to an inflation of type I error. Therefore, a Bonferroni correction is necessary. In our example the test results indicate that only the central branch is identified as saturated for category c4. The output we receive in the file \texttt{example\_aa.fasta\_c4\_0.01.satute.csv} is as follows:

        | edge                       | ... | **z_score** | ... | **z_alpha_bonferroni_corrected** | **decision_bonferroni_corrected** | ... | branch_length | number_of_sites | rate_category |
        |----------------------------|-----|-------------|-----|----------------------------------|---------------------------------------------|-----|---------------|-----------------|---------------|
        | ...                        | ... | ...         | ... | ...                              | ...                                         | ... | ...           | ...             |    c4         |
        | (Node5*, Node4*)           | ... | **12.9131** | ... | **3.5293**                       | **Informative**                             | ... | 0.2150928463  |  234            |    c4         |
        | (Node8*, Node4*)           | ... | **1.4587**  | ... | **3.6047**                       | **Saturated**                               | ... | 8.6358874072  |  234            |    c4         |
        | (Node6*, Node5*)           | ... | **13.9859** | ... | **3.384**                        | **Informative**                             | ... | 0.11624130057072  |  234            |    c4         |
        | ...                        | ... | ...         | ... | ...                              | ...                                         | ... | ...           | ...             |    c4         |
  
        Again, the central branch is saturated.

## Further Commands

Given a path to an IQ-Tree executable, SatuTe runs IQ-Tree with default options to generate the required data for the test of saturation, namely  a multiple sequence alignment (MSA), a model of sequence evolution with its parameters, and a phylogenetic tree. The following commands can be used to generate the same results using other SatuTe options.

Here, use the SatuTe mode MSA:

```bash
satute -msa ./data/example_aa.fasta -iqtree path_to_iqtree_exe -alpha 0.01
```
