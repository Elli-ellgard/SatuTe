import os
import shutil
import subprocess
import glob
from pathlib import Path


def print_info(suffix):
    print(" ============= MODI MSA ====================")
    print("")
    print(suffix)
    print("---------------------------------------------------------")


def create_dir_copy_with_suffix(source_dir, suffix) -> Path:
    suffix_dir = source_dir + "_" + suffix.replace(" ", "_").replace(":", "")
    if os.path.exists(suffix_dir):
        shutil.rmtree(suffix_dir)
    os.makedirs(suffix_dir, exist_ok=True)
    shutil.copytree(source_dir, suffix_dir, dirs_exist_ok=True)
    return suffix_dir


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
    dir_path = "./test/Clemens/toy_example_JC"
    msa = "toy_example_ntaxa_7_run_5-alignment.phy"
    path_iqtree = "iqtree2"
    path_python = "python3"

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
