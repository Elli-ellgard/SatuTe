# Satute
Satute is a tool to test for branch saturation in a phylogeny. As an intuition, a branch is saturated when it is too long, although sometimes the cause is more technical than just branch length.

Satute uses IQ-TREE for its core operations (http://www.iqtree.org) and provides a command-line interface for running saturation tests in a streamlined and manageable way. It also contains additional functionalities specific to saturation testing.

For the saturation test, we assume that the model of sequence evolution is reversible and stationary.

## Installation

Clone this repository to your local machine and navigate to the root directory. Satute is written in Python, so you need to have Python3 installed on your machine to use it. IQ-TREE should also be installed, and its executable should be within your system's PATH.

## Workflow

The minimum input required to run Satute is an alignment file in phylip or fasta format. IQ-TREE reconstructs an evolutionary tree and then Satute will tell us, for every branch of this tree, whether this branch is saturated or not. If the sites have different evolutionary speeds (e.g. due to the Gamma model), then Satute will separate the alignment sites into the rate categories that likely generated them.

The final output is a result file listing, for every branch, its saturation status. If there are different rate categories, then one result file is generated for each category.

The statistic we use to test for branch saturation is the coherence of the branch, which has a value typically between 0 (saturated) and 1 (unsaturated).

## Usage

Satute provides a command-line interface for easy usage. After installation, you can use the following command to run Satute:

```bash
python satute_cli.py -dir ./test/octo-kraken-msa-test -iqtree iqtree -model GTR+F
```
This command will start a saturation test on the data in the `test/octo-kraken-test` directory, using the IQ-TREE executable located at `./iqtree/bin/iqtree`.

or
```bash
python satute_cli.py -tree /path/to/mytree.nwk -msa /path/to/myalignment.fasta -iqtree /path/to/iqtree
```

This command will run the `satute.py` script with the following input:

- `-tree /path/to/mytree.nwk`: Path to the input tree file. Replace `/path/to/mytree.nwk` with the actual path to your Newick formatted tree file.
- `-msa /path/to/myalignment.fasta`: Path to the Multiple Sequence Alignment (MSA) file. Replace `/path/to/myalignment.fasta` with the actual path to your MSA file.
- `-iqtree /path/to/iqtree`: Path to the IQ-TREE executable file. Replace `/path/to/iqtree` with the actual path to your IQ-TREE executable.

The `-dir` flag can alternatively be used to specify a single directory where all data files are located, but it should not be used in combination with the `-tree` and `-msa` flags.


TODO: List all four possibilities (MSA, MSA+model, MSA +tree, MSA+tree+model)
## Parameters

Sure, here's how you could add the parameters to your markdown file:

## Satute Parameters

Here is the list of parameters you can use with Satute:

- `-dir <file_name>`: Path to the input directory containing the data files. Default is the current directory.

- `-tree <file_name>`: Path to the input tree file. Default is the input directory.

- `-msa <file_name>`: Path to the Multiple Sequence Alignment (MSA) file. Default is the input directory.

- `-iqtree <file_name>`: Path to the IQ-TREE executable file. The default path is assigned during the initialization of the class.

- `-model <str>`: Model name string (e.g. GTR+F+I+G). All reversible substitution models supported by IQ-TREE can be selected (see [IQ-TREE Substitution Models](http://www.iqtree.org/doc/Substitution-Models)). Default is assigned during the initialization of the class.

- `-nr <num>`: Number of rate categories for the Gamma model. Default is assigned during the initialization of the class.

- `-ufboot <num>`: Number of replicates for ultrafast bootstrap (must be >= 1000). Default is assigned during the initialization of the class.

- `-boot <num>`: Number of replicates for bootstrap + ML tree + consensus tree. Default is assigned during the initialization of the class.

You can specify these parameters when running the Satute command, with the values following the respective flags. Default values are used if not specified. 

## Methods

Satute provides the following methods:

- `parse_input()`: Parses command-line arguments.
- `write_log()`: Writes a log file that includes parameter values and other information.
- `run()`: Main entry point for running the Satute command-line tool.
- `run_iqtree_with_arguments(file_argument, extra_arguments=[])`: Runs IQ-TREE with given arguments and extra arguments.
- `check_input()`: Checks if the input directory exists and contains the required files.
- `find_file(suffixes)`: Finds file in input directory with given suffixes.
- `extract_best_model_from_log_file(file_path)`: Extracts best model from IQ-TREE log file.
- `file_exists(file_path)`: Checks if a file exists.

## Visualisation
To visualize the results of Statute, run the following bash script with an absolute path of the directory, where the results files of statute and the msa and tree file is present, as argument

```bash
cd visualisation
./script-annotate-trees-from-satute.sh $(pwd)/../test/visual_example/
```
The bash script executes the  R script
   1. script-annotate-tree-satute.R: Creation of newick string and nexus file with all meta data
   2. script-visualize-satute-tree.R: Generating three plots using delta value, p-value and branch status


## Satute Project To-Do List
1. **Diagonalisation in Saturation Test**
   - Remove the process of diagonalisation from the saturation test.

2. **Saturation Test Environment**
   - Move all pre-calculation processes out of the saturation test environment to streamline the workflow.

3. **IQ-TREE Execution for Each Clade**
   - Discuss the need for running IQ-TREE for every clade and the issue of every file being created twice. This could potentially lead to redundancy and could slow down the workflow.

4. **Dimensionality of Substitution Vectors**
   - Discuss the possibility that dimensions of substitution vectors can be more than 4. This could affect the calculation and interpretation of results.

5. **Calculation of R, State Frequencies and Q Tensor**
   - Move the calculation of R, state frequencies and Q tensor before the saturation test and perform it in memory. This change could speed up the overall processing time.

6. **Handling of Results File**
   - Rewrite the method of handling results files for more efficient and organized data management. Currently, the handling method may not be optimal and could be improved.

7. **Command Line Interface and Function Calls**
   - Currently, Satute is designed to be used as a command-line tool which can make calls from another Python script challenging. Consider implementing a way for Satute to be easily called from another Python script.
    
8. **Calculation of test according to the non-dominant, non-zero eigenvalues**
9. **Integrate other arguments like seed and significance level, remove c_s value**
10. **Check calculation of the delta value, especially for multiplicity >1**

By addressing these tasks, we can improve the efficiency and flexibility of the Satute project, allowing it to be utilized more effectively in a variety of workflows. Remember that productive collaboration and regular code reviews will significantly contribute to achieving these goals.

## Contributions

Contributions to Satute are welcome! Please fork this repository and submit a pull request with your changes. If you find any issues, feel free to report them on the issue tracker.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
