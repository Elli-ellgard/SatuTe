# Satute
Satute is a asymptotic test for branch saturation tool that utilizes IQ-TREE for its core operations. It provides a command-line interface for running saturation tests in a more streamlined and manageable way. While it heavily interacts with IQ-TREE, it is not merely a wrapper but incorporates additional functionalities specific to saturation testing. For the saturation test, we assume as model of sequence evolution a general reversible model in equilibrium.

## Installation

Clone this repository to your local machine and navigate to the root directory. Satute is written in Python, so you need to have Python installed on your machine to use it. IQ-TREE should also be installed, and its executable should be within your system's PATH.

## Usage

Satute provides a command-line interface for easy usage. After installation, you can use the following command to run Satute:

```bash
python satute_cli.py -dir ./test/octo-kraken-msa-test -iqtree iqtree -model GTR+F
```

This command will start a saturation test on the data in the `test/octo-kraken-test` directory, using the IQ-TREE executable located at `./iqtree/bin/iqtree`.

TODO: List all four possibilities (MSA, MSA+model, MSA +tree, MSA+tree+model)

## Parameters

- `-iqtree` (str): Path to the IQ-TREE executable. Default is `"iqtree"`.
- `-dir` (str): Path to the input directory containing the data files. Default is `"./"`. The directory should contain at least a PHYLIP/FASTA alignment file.
- `-m` (str): Model name string (e.g. GTR+F+I+G). Default is `None`.  All reversible substitution models supported by IQ-TREE can be selected (see http://www.iqtree.org/doc/Substitution-Models)
- `-nr` (int): Number of rate categories for the Gamma model. Default is `None`.
- `-o` (str): Prefix for the output files. Default is `None`.
- `-ufboot` (int): Number of replicates for ultrafast bootstrap (must be >= 1000). Default is `None`.
- `-boot` (int): Number of replicates for bootstrap + ML tree + consensus tree. Default is `None`.

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

## Contributions

Contributions to Satute are welcome! Please fork this repository and submit a pull request with your changes. If you find any issues, feel free to report them on the issue tracker.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
