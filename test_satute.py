import os
import subprocess

iqtree = "iqtree2"


def execute_command(command):
    subprocess.run(command, shell=True)


def process_directory(directory):
    # execute_command("source ./env/bin/activate")
    execute_command(
        f"python3 satute_cli.py -iqtree {iqtree} -msa {directory} -alpha 0.05"
    )


# Test cases from Clemens
test_cases = [
    # ("./Clemens/example_1", "ENSG00000087460_GNAS.fasta"),
    # (
    #     "./Clemens/example_2",
    #     "sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy",
    # ),
    (
        "./Clemens/example_3",
        "sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    ),
    # ("./Clemens/example_4", "example.txt"),
    # ("./Clemens/example_sym_1", "ENSG00000119574_ZBTB45.fasta"),
    # ("./Clemens/example_sym_2", "ENSG00000138316_ADAMTS14.fasta"),
    # ("./Clemens/toy_example_GTR+G4", "toy_example_ntaxa_7_run_1-alignment.phy"),
    # ("./Clemens/toy_example_JC", "toy_example_ntaxa_7_run_5-alignment.phy"),
]

for dir_path, msa in test_cases:
    print(f"Processing {dir_path}/{msa}")
    process_directory(f"{dir_path}/{msa}")
