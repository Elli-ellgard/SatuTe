import os
import subprocess
import time


import sys
# Add the directory containing the satute module to the Python path
sys.path.append('/home/elgert/Desktop/Cassius/version_2024_08_19')

import satute.cli

def find_file_with_suffix_in_directory(suffix, source_folder):
    # Iterate over files in the source folder
    for filename in os.listdir(source_folder):
        if filename.endswith(suffix):
            return os.path.join(source_folder, filename)
    return None

def run_external_command(command_args):
    subprocess.run(command_args)

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
    if '-quiet' not in args:
        args.append('-quiet')

    # Add a small delay to ensure files are written
    time.sleep(2)

    # Call the main function of satute.cli with the arguments
    return satute.cli.main(args)

def  run_satute_for_directory(folder_path, alpha):

    run_satute(
        [
            "-dir",
            folder_path,
            "-alpha",
            alpha,
        ]
    )
    #run_external_command(arguments)
        

def clean_paths_satute(input_dir, path_prefix):
    def clean_line(line, base_dir):
        # Function to clean a single line by removing the path prefix
        return line.replace(base_dir, '')


    # Find the .satute and .satute.log files
    satute_file = find_file_with_suffix_in_directory(".satute", input_dir)
    log_file = find_file_with_suffix_in_directory(".satute.log", input_dir)

    # Process the .satute file
    if satute_file:
        with open(satute_file, 'r') as file:
            lines = file.readlines()
        with open(satute_file, 'w') as file:
            for line in lines:
                file.write(clean_line(line, path_prefix))
        print(f"Processed: {satute_file}")

    # Process the .satute.log file
    if log_file:
        with open(log_file, 'r') as file:
            lines = file.readlines()
        with open(log_file, 'w') as file:
            for line in lines:
                file.write(clean_line(line, path_prefix))
        print(f"Processed: {log_file}")



if __name__ == "__main__":

    # Set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"
  
    """ Set directories and input files """
    # Get the current working directory
    current_directory = os.getcwd()

    # Examples DNA alignments
    example1 = os.path.join(current_directory, "examples_dna/dir_true_tree")
    alpha = str(0.05)
    run_satute_for_directory(example1, alpha)
    clean_paths_satute(example1, current_directory)

    try:
        subprocess.run('mv *satute* ./SatuTe_results/', shell=True, check=True, cwd=example1)
        print("Files moved successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
    
    


    example2 = os.path.join(current_directory, "examples_dna/dir_ML_tree")
    alpha = str(0.01)
    run_satute_for_directory(example2, alpha)
    clean_paths_satute(example2, current_directory)
    try:
        subprocess.run('mv *satute* ./SatuTe_results/', shell=True, check=True, cwd=example2)
        print("Files moved successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")


    example3 = os.path.join(current_directory, "examples_aa/dir_ML_tree")
    alpha = str(0.01)
    run_satute_for_directory(example3, alpha)
    clean_paths_satute(example3, current_directory)
    try:
        subprocess.run('mv *satute* ./SatuTe_results/', shell=True, check=True, cwd=example3)
        print("Files moved successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

    #clean_paths_satute(f"{example1}/SatuTe_results/", current_directory)