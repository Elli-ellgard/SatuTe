import os
import shutil
import subprocess
import glob
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
    return suffix

def check_iqtree_files_exist(data_name, dest_dir_path, iqtree_options):
    iqtree_files_endings = [".bionj",".ckp.gz",".iqtree",".log",  ".mldist", ".model.gz", ".treefile"]
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


def test_1a(dir_path, msa, iqtree, python):
    suffix = "TEST 1a: only msa iqtree exists"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)

    run_external_command(
        [
            python,
            "satute_cli.py",
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(new_dir_path, msa),
            "-alpha",
            "0.05",
            "-asr",
        ]
    )


def test_1b(dir_path, msa, iqtree, python):
    suffix = "TEST 1b: only msa ufboot"
    print(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)

    run_external_command(
        [
            python,
            "satute_cli.py",
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(new_dir_path, msa),
            "-ufboot",
            "1000",
            "-alpha",
            "0.05",
            "-asr",
        ]
    )


def test_1c(dir_path, msa, iqtree, python):
    suffix = "TEST 1c: only msa boot"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            python,
            "satute_cli.py",
            "-iqtree",
            iqtree,
            "-msa",
            os.path.join(new_dir_path, msa),
            "-boot",
            "100",
            "-alpha",
            "0.05",
            "-asr",
        ]
    )


# New test functions
def test_1d(dir_path, msa, path_iqtree):
    suffix = "TEST 1d: only msa iqtree exists path to executable"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-iqtree",
            path_iqtree,
            "-msa",
            os.path.join(new_dir_path, msa),
            "-alpha",
            "0.05",
        ]
    )


def test_2a(path_iqtree):
    suffix = "TEST 2a: no msa  no dir => error"
    print_info(suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-model",
            "JC",
            "-alpha",
            "0.05",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_2b(dir_path, msa):
    suffix = "TEST 2b: only msa iqtree not exists => error"
    print_info(suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-iqtree",
            "iqtree3",
            "-msa",
            os.path.join(dir_path, msa),
            "-alpha",
            "0.05",
            "-asr",
        ]
    )


def test_3a(dir_path, msa, path_iqtree):
    suffix = "TEST 3a: msa + model number of rate categories set"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC+G4",
            "-alpha",
            "0.05",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_3c(dir_path, msa, path_iqtree):
    suffix = "TEST 3c: msa + model with parameter + shape parameter of Gamma model"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "TIM2{4.39/5.30/12.1}+G{0.7}",
            "-alpha",
            "0.05",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_3c2(dir_path, msa, path_iqtree):
    suffix = "TEST 3c: msa + model with parameter"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "TIM2{4.39/5.30/12.1}+R2{0.7/0.1/0.3/0.9}",
            "-alpha",
            "0.05",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_3d(dir_path, msa, path_iqtree):
    suffix = "TEST 3d: msa + model option ufboot"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC+G4",
            "-ufboot",
            "1000",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_3e(dir_path, msa, path_iqtree):
    suffix = "TEST 3e: msa + model option boot"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC+G4",
            "-boot",
            "100",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_4(dir_path, msa, path_iqtree):
    suffix = "TEST 4: msa + model + nr"
    print_info(suffix)
    new_dir_path = create_dir_copy_with_suffix(dir_path, suffix)
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC+G4",
            "-category",
            "2",
            "-alpha",
            "0.05",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_5(dir_path, msa, iqtree):
    print_info("TEST 5: msa + model + tree")
    new_dir_path = create_dir_copy_with_suffix(dir_path, "TEST 5: msa + model + tree")

    # Clean and prepare directory
    clean_and_prepare_dir(new_dir_path, msa)

    # First part of the test - creating output
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-iqtree",
            iqtree,
        ]
    )

    # Clean and prepare directory again
    clean_and_prepare_dir(new_dir_path, msa)

    # Second part of the test
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC+G4",
            "-tree",
            os.path.join(new_dir_path, msa + ".treefile"),
            "-iqtree",
            iqtree,
        ]
    )


def test_5_ufboot_error(dir_path, msa, iqtree):
    print_info("TEST 5: msa + model + tree, option ufboot => error")

    # Clean and prepare directory for the test
    clean_and_prepare_dir(dir_path, msa)

    # Create initial output
    run_external_command(
        ["python", "satute_cli.py", "-msa", os.path.join(dir_path, msa)]
    )

    # Clean and prepare directory again
    clean_and_prepare_dir(dir_path, msa)

    # Run the test command expected to produce an error
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(dir_path, msa),
            "-model",
            "JC+G4",
            "-tree",
            os.path.join(dir_path, msa + ".treefile"),
            "-ufboot",
            "1000",
        ]
    )


def test_6_directory_error(dir_path):
    print_info("TEST 6: directory does not exist => error")

    # Construct the non-existing directory path
    non_existing_dir_path = dir_path + "_bla_bla_bla"

    # Run the command expected to produce an error
    run_external_command(["python", "satute_cli.py", "-dir", non_existing_dir_path])


def test_7_dir_missing_siteprob(original_dir_path, msa, iqtree_path):
    print_info("TEST 7: dir missing siteprob")

    # Create a new directory with a suffix
    new_dir_path = create_dir_copy_with_suffix(
        original_dir_path, "TEST_7_dir_missing_siteprob"
    )

    # Run IQ-TREE with the specified model in the new directory
    run_external_command(
        [iqtree_path, "-s", os.path.join(new_dir_path, msa), "-m", "JC+R2", "-quiet"]
    )

    # Run satute_cli.py with the new directory option
    run_external_command(["python", "satute_cli.py", "-dir", new_dir_path])


def test_8_dir_msa(dir_path, msa):
    print_info("TEST 8: dir + msa")

    # Clean and prepare directory for the test
    clean_and_prepare_dir(dir_path, msa)

    # Run the first part of the test
    run_external_command(
        ["python", "satute_cli.py", "-msa", os.path.join(dir_path, msa)]
    )

    # Remove any satute tree files generated in the directory
    for file in glob.glob(os.path.join(dir_path, f"{msa}*satute*")):
        os.remove(file)

    # Run the second part of the test
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-dir",
            dir_path,
            "-msa",
            os.path.join(dir_path, msa),
        ]
    )


def test_5_edge_option(dir_path, msa, path_iqtree):
    print_info("TEST 5: msa + model + tree option edge")

    # Create a new directory with a suffix and run the first part of the test
    new_dir_path = create_dir_copy_with_suffix(dir_path, "TEST_5_edge_option")

    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-iqtree",
            path_iqtree,
        ]
    )

    # Run the second part of the test with model, tree and edge option
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_dir_path, msa),
            "-model",
            "JC",
            "-tree",
            os.path.join(new_dir_path, msa + ".treefile"),
            "-edge",
            "(Node2*, Node1*)",
            "-iqtree",
            path_iqtree,
        ]
    )


def test_5_dir_option_edge(dir_path, msa):
    print_info("TEST 5: dir option edge")

    # Create a new directory with a suffix and run the first part of the test
    new_dir_path = create_dir_copy_with_suffix(dir_path, "TEST_5_dir_option_edge")
    run_external_command(
        ["python", "satute_cli.py", "-msa", os.path.join(new_dir_path, msa)]
    )

    # Run the second part of the test with dir and edge options
    run_external_command(
        ["python", "satute_cli.py", "-dir", new_dir_path, "-edge", "(Node2*, Node1*)"]
    )


def test_10a_invalid_category(dir_path, msa, iqtree):
    print_info("TEST 10a: msa + model, option invalid category => error")

    # Running the command with invalid category option
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(dir_path, msa),
            "-model",
            "JC+G4",
            "-category",
            "5",
            "-iqtree",
            iqtree,
            "-verbose",
        ]
    )


def test_10c_valid_category(dir_path, msa, iqtree):
    suffix = "TEST 10c: directory option valid category"
    print_info(suffix)
    # Remove any satute tree files generated in the directory
    new_directory = create_dir_copy_with_suffix(
        dir_path,
        "TEST 10c: directory option valid category".replace(" ", "_")
        .replace(":", "")
        .lower(),
    )

    # Run the first part of the test to create output
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-msa",
            os.path.join(new_directory, msa),
            "-model",
            "GTR+G6",
            "-iqtree",
            iqtree,
            "-asr"            
        ]
    )

    # Run the second part of the test with a valid category option
    run_external_command(
        [
            "python",
            "satute_cli.py",
            "-dir",
            new_directory,
            "-category",
            "1",
        ]
    )


def main():

    # set paths to IQ-TREE and Python executable
    path_iqtree = "iqtree2"
    path_python = "python3"



    dir_path = "./test/Clemens/toy_example_JC"
    msa = "toy_example_ntaxa_7_run_5-alignment.phy"


    test_1a(dir_path, msa, path_iqtree, path_python)
    test_1b(dir_path, msa, path_iqtree, path_python)
    test_1c(dir_path, msa, path_iqtree, path_python)

    # test_2a(path_iqtree, path_python)
    # test_2b(dir_path, msa, path_python)

    # # Run new tests
    # test_3c(dir_path, msa, path_iqtree, path_python)
    # test_3c2(dir_path, msa, path_iqtree, path_python)
    # test_3d(dir_path, msa, path_iqtree, path_python)
    # test_3e(dir_path, msa, path_iqtree, path_python)

    # test_4(dir_path, msa, path_iqtree, path_python)

    # # Run new test
    # test_5(dir_path, msa, path_iqtree, path_python)
    # test_5_edge_option(dir_path, msa, path_iqtree, path_python)

    # test_6_directory_error(dir_path, path_python)

    # test_7_dir_missing_siteprob(dir_path, msa, path_iqtree, path_python)

    # test_8_dir_msa(dir_path, msa, path_python)

    # test_10a_invalid_category(dir_path, msa, path_iqtree, path_python)
    # test_10c_valid_category(dir_path, msa, path_iqtree, path_python)


if __name__ == "__main__":
    main()
