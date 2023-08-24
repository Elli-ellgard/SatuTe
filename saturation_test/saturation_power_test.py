import subprocess
import random
import string
import os


iqtree = "iqtree"
from ete3 import Tree


def fetch_files_with_prefix(directory, prefix):
    """
    Fetch all files in a directory that start with a specific prefix.

    Parameters:
    - directory: The path to the directory to search.
    - prefix: The file prefix to match.

    Returns:
    - A list of filenames that match the prefix.
    """

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files based on the prefix
    matching_files = [f for f in all_files if f.startswith(prefix)]

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


import os


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


def process(model, sequence_length, seed, alignment_num):
    prefix = f"{model}-{sequence_length}"
    # execute_command("source ./env/bin/activate")
    generate_tree(10)
    execute_command(
        f"iqtree --alisim {model}-{sequence_length} -t random_generated.tree --length {sequence_length} --seqtype DNA --num-alignments {alignment_num} -m {model}"
    )

    all_files = fetch_files_with_prefix(".", f"{model}-{sequence_length}")
    # Filter files based on the prefix
    files = [f for f in all_files if f.startswith(prefix)]    
    
    execute_command(
         f"python3 ./../satute_cli.py -iqtree {iqtree} -msa JC-1000_1.phy -tree ./random_generated.tree -model {model} -alpha 0.05"
    )
    
    # execute_command(
    #      f"python3 ./../satute_cli.py -iqtree {iqtree} -msa JC-1000_2.phy -tree random_generated.tree -model {model} -alpha 0.05"
    # )
    # Example usage:
    

if __name__ == "__main__":
    # model_set = ["GTR", "GTR+G4", "GTR+G8", "JC", "JC+G4", "JC+G8"]
    prefix = "JC-1000"
    delete_files_with_prefix("./", prefix)
    process("JC", 100, 1, 3)
   