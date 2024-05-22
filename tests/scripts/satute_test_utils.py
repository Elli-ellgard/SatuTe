import os
import shutil
import subprocess
from pathlib import Path

def print_test_name(suffix):
    print("")
    print("---------------------------------------------------------")
    print(suffix)
    print("---------------------------------------------------------")

def print_colored_message(message, color_code):
    print(f"\033[{color_code}m{message}\033[0m")

def create_destination_dir(source_dir, suffix, base_dir="../test_results/") -> Path:
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
    
    missing_files = [file for file in file_list if not (directory_path / file).is_file()]
    
    if missing_files:
        print(f"The following {name} files are missing in the directory '{directory}': {', '.join(missing_files)}")
        return False
    else:
        print(f"All {name} files are present in the directory '{directory}'.")
        return True


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
    return suffix

def check_iqtree_files_exist(data_name, dest_dir_path, iqtree_options):
    iqtree_files_endings = [".bionj",".ckp.gz",".iqtree",".log",  ".mldist" , ".treefile"]
    if len(iqtree_options):
        iqtree_files_endings.extend(additional_iqtree_files(iqtree_options))
    files_to_check = [ data_name + suffix for suffix in iqtree_files_endings]
    return check_files_exist(dest_dir_path, files_to_check, "IQ-TREE")


def check_satute_files(data_name, dest_dir_path, categories, alpha, asr):
    file_endings = [f"_{alpha}.satute.log"]
    suffix = [".satute.csv",".satute.nex"]
    if asr: 
        suffix.append(".satute.asr.csv")
    if len(categories):
        for rate in categories:
            category_endings = [f"_p{rate}_{alpha}{x}" for x in suffix]
            file_endings.extend(category_endings)
    else:
        file_endings.extend([f"_single_rate_{alpha}{x}" for x in suffix])
        
    files_to_check = [ data_name + suffix for suffix in file_endings]  
    return check_files_exist(dest_dir_path, files_to_check, "Satute")



def run_external_command(command_args):
    subprocess.run(command_args)


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




