import os
import shutil
import subprocess
import time
from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree
import satute.cli


def filter_files_by_last_name(filenames, suffix):
    """
    Filters the filenames that end with any of the specified suffixes.

    Args:
    filenames (list): List of filenames to filter.
    suffixes (list): List of suffixes to match filenames against.

    Returns:
    list: List of filenames that match the suffixes.
    """
    filtered_list = []
    for file in filenames:
        splitted_file_name = file.split(".")

        if splitted_file_name[-1] != suffix:
            filtered_list.append(file)

    return filtered_list


def create_destination_dir(
    source_dir, suffix, base_dir="./tests/test_results/"
) -> Path:
    source_path = Path(source_dir)
    base_path = Path(base_dir)

    # Check if the source directory exists
    if not source_path.is_dir():
        raise ValueError("Source directory does not exist.")

    # Create a suffix for the new directory name
    suffix_cleaned = suffix.replace(" ", "_").replace(":", "")

    # Create the destination directory path
    dest_dir = base_path / (source_path.name + "_" + suffix_cleaned)

    # Remove destination directory if it already exists
    if dest_dir.is_dir():
        shutil.rmtree(dest_dir)

    # Create the destination directory
    dest_dir.mkdir(parents=True, exist_ok=True)

    return dest_dir


def copy_files_to_dest_dir(source_dir, dest_dir, files_to_copy):
    source_path = Path(source_dir)
    dest_path = Path(dest_dir)

    # Check if both source and destination directories exist
    if not source_path.is_dir() or not dest_path.is_dir():
        raise ValueError("Source or destination directory does not exist.")

    # Copy only the specified files from source to destination
    for file_name in files_to_copy:
        source_file_path = source_path / file_name
        dest_file_path = dest_path / file_name

        # Check if the file exists in the source directory
        if source_file_path.is_file():
            shutil.copy2(source_file_path, dest_file_path)
        else:
            print(f"Warning: File '{file_name}' not found in the source directory.")


def check_files_exist(directory, file_list, name):
    directory_path = Path(directory)

    # Check if the directory exists
    if not directory_path.is_dir():
        raise ValueError(f"Directory '{directory}' does not exist.")

    missing_files = [
        file for file in file_list if not (directory_path / file).is_file()
    ]

    if missing_files:
        print(
            f"The following {name} files are missing in the directory '{directory}': {', '.join(missing_files)}"
        )
        return False
    else:
        print(f"All {name} files are present in the directory '{directory}'.")
        return True


def check_files_exist_dir(directory_path, filenames, description):
    """
    Checks if the specified files exist in the directory.

    Args:
    directory_path (str): Path to the directory.
    filenames (list): List of filenames to check.
    description (str): Description of the file type (for reporting).

    Returns:
    dict: Dictionary with filenames as keys and existence as boolean values.
    """
    file_existence = {}
    for filename in filenames:
        file_path = os.path.join(directory_path, filename)
        file_existence[filename] = os.path.isfile(file_path)
    return file_existence


def additional_iqtree_files(options):
    suffix = []
    if "ufboot" in options:
        suffix.extend([".contree", ".splits.nex"])
    elif "boot" in options:
        suffix.extend([".boottrees", ".contree"])
    elif "wspr" in options:
        suffix.append(".siteprob")
    elif "m" in options:
        suffix.append(".model.gz")
    elif "alnifo" in options:
        suffix.append(".alninfo")
    return suffix


def check_iqtree_files_exist(data_name, dest_dir_path, iqtree_options):
    iqtree_files_endings = [
        ".bionj",
        ".ckp.gz",
        ".iqtree",
        ".log",
        ".mldist",
        ".treefile",
    ]
    if len(iqtree_options):
        iqtree_files_endings.extend(additional_iqtree_files(iqtree_options))
    files_to_check = [data_name + suffix for suffix in iqtree_files_endings]
    return check_files_exist(dest_dir_path, files_to_check, "IQ-TREE")


def check_satute_files(data_name, dest_dir_path, categories, alpha, asr):
    file_endings = [f"_{alpha}.satute.log"]
    suffix = [".satute.csv", ".satute.nex"]
    if asr:
        suffix.append(".satute.asr.csv")
    if len(categories):
        for rate in categories:
            category_endings = [f"_c{rate}_{alpha}{x}" for x in suffix]
            file_endings.extend(category_endings)
    else:
        file_endings.extend([f"_single_rate_{alpha}{x}" for x in suffix])

    files_to_check = [data_name + suffix for suffix in file_endings]
    return check_files_exist(dest_dir_path, files_to_check, "Satute")


def check_satute_files_dir(dest_dir_path, categories, alpha, asr):
    """
    Constructs and checks the existence of expected files with given suffixes for 'satute'.

    Args:
    dest_dir_path (str): Path to the destination directory where files are expected.
    categories (list): List of category rates.
    alpha (str): A suffix to append to the file names.
    asr (bool): Whether to include the ASR suffix.

    Returns:
    dict: Dictionary with filenames as keys and existence as boolean values.
    """
    file_endings = [f"_{alpha}.satute.log"]
    suffixes = [".satute.csv", ".satute.nex"]

    if asr:
        suffixes.append(".satute.asr.csv")

    if categories:
        for rate in categories:
            category_endings = [f"_c{rate}_{alpha}{suffix}" for suffix in suffixes]
            file_endings.extend(category_endings)
    else:
        file_endings.extend([f"_single_rate_{alpha}{suffix}" for suffix in suffixes])

    # Construct the filenames based on the suffixes
    files_to_check = file_endings

    # Check for file existence
    return check_files_exist_dir(dest_dir_path, files_to_check, "Satute")


def run_external_command(command_args):
    subprocess.run(command_args, check=True)


def run_satute(args):
    """
    Runs the satute CLI with the provided arguments, ensuring the 'quiet'
    argument is set, and adds a small delay to ensure files are written.

    Args:
        args (list or None): Command-line arguments to pass to satute.cli.main. If None, an empty list is used.

    Returns:
        bool: The result of satute.cli.main(args).
    """
    if args is None:
        args = []

    if not isinstance(args, list):
        raise TypeError("args must be a list or None")

    # Ensure '--quiet' is in the arguments
    if "-quiet" not in args:
        args.append("-quiet")

    # Add a small delay to ensure files are written
    time.sleep(2)

    # Call the main function of satute.cli with the arguments
    return satute.cli.main(args)


def clean_and_prepare_dir(dir_path, msa):
    pdir = os.path.dirname(dir_path)
    msa_path = os.path.join(dir_path, msa)
    iqtree_path = msa_path + ".iqtree"

    if os.path.exists(iqtree_path):
        shutil.move(msa_path, pdir)
        if os.path.exists(os.path.join(dir_path, msa + ".treefile")):
            shutil.move(os.path.join(dir_path, msa + ".treefile"), pdir)
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path)
        shutil.move(os.path.join(pdir, msa), dir_path)
        if os.path.exists(os.path.join(pdir, msa + ".treefile")):
            shutil.move(os.path.join(pdir, msa + ".treefile"), dir_path)


def list_filenames_in_directory(directory_path):
    """
    Lists all filenames in the given directory and its subdirectories.

    Args:
    directory_path (str): Path to the directory.

    Returns:
    list: List of filenames.
    """
    filenames = []

    for root, dirs, files in os.walk(directory_path):
        for file in files:
            filenames.append(file)

    return filenames


def delete_taxon_from_msa(
    input_file: str, output_file: str, taxon_id: str, format: str = "phylip-sequential"
):
    """
    Delete a taxon from a multiple sequence alignment.

    Args:
        input_file (str): Path to the input MSA file.
        output_file (str): Path to the output MSA file with the taxon removed.
        taxon_id (str): The identifier of the taxon to remove.
        format (str): Format of the MSA file (default: "fasta").
    """
    # Read the MSA from the input file
    alignment = AlignIO.read(input_file, format)

    # Filter out the sequences with the specified taxon_id
    filtered_alignment = MultipleSeqAlignment(
        record for record in alignment if record.id != taxon_id
    )

    # Write the filtered alignment to the output file
    AlignIO.write(filtered_alignment, output_file, format)

    print(
        f"Taxon '{taxon_id}' has been removed. The new MSA is saved to '{output_file}'."
    )


def delete_taxon_from_tree(newick_file: str, taxon: str):
    """
    Delete a taxon from a phylogenetic tree and overwrite the Newick file with the updated tree.

    Args:
        newick_file (str): Path to the Newick file containing the tree.
        taxon (str): The name of the taxon to delete.

    Raises:
        FileNotFoundError: If the Newick file does not exist.
        ValueError: If the taxon is not found in the tree.
    """
    # Load the tree from the Newick file
    try:
        tree = Tree(newick_file)
    except FileNotFoundError as e:
        print(f"File not found: {newick_file}")
        raise e

    # Find and delete the taxon
    node = tree.search_nodes(name=taxon)
    if not node:
        raise ValueError(f"Taxon '{taxon}' not found in the tree.")

    node[0].delete()

    # Write the updated tree back to the Newick file
    with open(newick_file, "w") as f:
        f.write(tree.write(format=1))  # format=1 for Newick with branch lengths

    print(f"Taxon '{taxon}' has been removed and the tree has been updated.")


def find_file_with_suffix(gene_name, suffix, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if gene_name in filename and filename.endswith(suffix):
            return os.path.join(source_folder, filename)
    return None
