import subprocess
import random
import string
import os
import shutil

from ete3 import Tree


iqtree = "iqtree"


def fetch_files_with_prefix_and_extension(directory, prefix, extension):
    """
    Fetch all files in a directory that start with a specific prefix and have a specific file type (extension).

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The file prefix to match.
    - extension: The file extension to match (e.g., ".txt" for text files).

    Returns:
    - A list of filenames that match the prefix and extension.
    """

    # Ensure the extension starts with a dot
    if not extension.startswith("."):
        extension = "." + extension

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files based on the prefix and extension
    matching_files = [
        f for f in all_files if f.startswith(prefix) and f.endswith(extension)
    ]

    return matching_files


def execute_command(command):
    subprocess.run(command, shell=True)


def generate_tree(n, file_name="random_generated.tree"):
    cell_names = [f"{random.choice(string.ascii_uppercase)}_cell-{i}" for i in range(n)]

    t = Tree()

    t.populate(
        size=n,
        names_library=cell_names,
        reuse_names=False,
        branch_range=(0.1, 1.0),
        support_range=(0.5, 1.0),
    )

    file_name = "random_generated.tree"

    with open(f"{file_name}", "w") as f:
        f.write(t.write(format=1))


def delete_files_with_prefix(directory, prefix):
    """
    Delete all files in a directory that start with a specific prefix.

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The file prefix to match.
    """
    # List all files in the directory
    all_files = os.listdir(directory)
    # Filter files based on the prefix and delete them
    for f in all_files:
        if f.startswith(prefix):
            file_path = os.path.join(directory, f)
            os.remove(file_path)
            print(f"Deleted: {file_path}")


def move_files_to_directory(prefix, alignment_num):
    # Assuming that the alignments and trees have specific extensions, e.g., ".phy" and ".tree"
    # Adjust these as needed
    alignment_extension = ".phy"

    for i in range(1, alignment_num + 1):
        alignment_file = f"{prefix}_{i}{alignment_extension}"

        # Create a new directory for the alignment and tree
        new_dir = f"{prefix}_{i}"
        os.makedirs(new_dir, exist_ok=True)

        if os.path.exists(os.path.join(new_dir, os.path.basename(alignment_file))):
            os.remove(os.path.join(new_dir, os.path.basename(alignment_file)))

        # Then line 96 remains as:
        shutil.move(alignment_file, new_dir)

    print(f"Moved {alignment_num} alignments and trees to individual directories.")




def execute_satute_cli_for_experiment_folders(main_directory, iqtree_executable):
    """
    Execute the satute_cli.py script for each experiment subdirectory in the main directory.

    Parameters:
    - main_directory (str): Path to the main directory containing experiment subdirectories.
    - iqtree_executable (str): Path to the IQ-TREE executable.
    - model (str): Substitution model to be used.

    Returns:
    - None
    """

    # List all subdirectories in the main directory
    subdirectories = [
        os.path.join(main_directory, d)
        for d in os.listdir(main_directory)
        if os.path.isdir(os.path.join(main_directory, d))
    ]

    for subdir in subdirectories:
        # Find the alignment file in the subdirectory (assuming they have .fasta extension)
        alignment_files = [f for f in os.listdir(subdir) if f.endswith(".phy")]

        # If there's an alignment file, process it
        if alignment_files:
            alignment_path = os.path.join(subdir, alignment_files[0])

            # Construct the command
            cmd = [
                "python3",
                "./../satute_cli.py",
                "-iqtree",
                iqtree_executable,
                "-msa",
                alignment_path,
            ]

            # Execute the command
            subprocess.run(cmd)


# Updating the process function to include the move step
def power_test_one_process(model, sequence_length, alignment_num):
    prefix = f"{model}-{sequence_length}"
    # execute_command("source ./env/bin/activate")
    generate_tree(10)
    os.chdir("./power_test_one")
    execute_command(
        f"iqtree --alisim {model}-{sequence_length} -t random_generated.tree --length {sequence_length} --seqtype DNA --num-alignments {alignment_num} -m {model}"
    )
    # Move the generated files to individual directories
    move_files_to_directory(prefix, alignment_num)
    execute_satute_cli_for_experiment_folders("./", iqtree)
    os.chdir("./..")


def power_test_two_process(model, alignment_num):
   
    os.chdir("./power_test_two")
    generate_tree(10, file_name="./random_generated.tree")

    for sequence_length in range(100, 1000, 100):
        prefix = f"{model}-{sequence_length}"
        execute_command(
            f"iqtree --alisim {prefix} -t random_generated.tree --length {sequence_length} --seqtype DNA --num-alignments {alignment_num} -m {model}"
        )
        move_files_to_directory(prefix, alignment_num)
        execute_satute_cli_for_experiment_folders("./", iqtree)


if __name__ == "__main__":
    # model_set = ["GTR", "GTR+G4", "GTR+G8", "JC", "JC+G4", "JC+G8"]
    prefix = "JC-100"
    # delete_files_with_prefix("./saturation_test", prefix)
    # power_test_one_process("JC", 100, 3)
    power_test_two_process("JC", 3)
