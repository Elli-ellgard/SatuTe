import subprocess

iqtree = "iqtree"


def execute_command(command):
    subprocess.run(command, shell=True)


def process_with_msa_tree_model(directory):
    # execute_command("source ./env/bin/activate")
    execute_command(
        f"python3 satute_cli.py -iqtree {iqtree} -msa {directory['msa_path']} -tree {directory['tree_file']} -model {directory['model']}"
    )


def process_without_model(directory):
    # execute_command("source ./env/bin/activate")
    execute_command(
        f"python satute_cli.py -iqtree {iqtree} -msa {directory['msa_path']} -tree {directory['tree_file']}"
    )


def process_without_model(directory):
    # execute_command("source ./env/bin/activate")
    execute_command(f"python3 satute_cli.py -iqtree {iqtree} -msa {directory}")


# Test cases from Clemens
test_cases_msa_tree_model = [
    # ("./Clemens/example_1", "ENSG00000087460_GNAS.fasta"),
    # (
    #     "./Clemens/example_2",
    #     "sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy",
    # ),
    # (
    #    "./Clemens/example_3",
    #    "sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    # ),
    # ("./Clemens/example_4", "example.txt"),
    # ("./Clemens/example_sym_2", "ENSG00000138316_ADAMTS14.fasta"),
    # ("./Clemens/toy_example_GTR+G4", "toy_example_ntaxa_7_run_1-alignment.phy"),
    {
        "msa_path": "./test/toy_example_JC/msa_files/toy_example_ntaxa_7_run_5-alignment.phy",
        "tree_file": "./test/toy_example_JC/tree_files/toy_example_ntaxa_7_run_5-alignment.phy.treefile",
        "model": "JC+G4",
    },
]

# Test cases from Clemens
test_cases_msa_tree_model = [
    # ("./Clemens/example_1", "ENSG00000087460_GNAS.fasta"),
    # (
    #     "./Clemens/example_2",
    #     "sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy",
    # ),
    # (
    #    "./Clemens/example_3",
    #    "sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
    # ),
    # ("./Clemens/example_4", "example.txt"),
    # ("./Clemens/example_sym_2", "ENSG00000138316_ADAMTS14.fasta"),
    # ("./Clemens/toy_example_GTR+G4", "toy_example_ntaxa_7_run_1-alignment.phy"),
    {
        "msa_path": "./test/toy_example_JC/msa_files/toy_example_ntaxa_7_run_5-alignment.phy",
        "tree_file": "./test/toy_example_JC/tree_files/toy_example_ntaxa_7_run_5-alignment.phy.treefile",
        "model": "JC+G4",
    },
]

test_case_without_model = [
    {
        "msa_path": "./test/toy_example_JC/msa_files/toy_example_ntaxa_7_run_5-alignment.phy",
        "tree_file": "./test/toy_example_JC/tree_files/toy_example_ntaxa_7_run_5-alignment.phy.treefile",
    },
]

# Test cases from Clemens
#test_cases_msa = [
#     "./Clemens/example_1/ENSG00000087460_GNAS.fasta",
#     "./Clemens/example_2/sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy",
#     "./Clemens/example_3/sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta",
#     "./Clemens/example_4/example.txt",
#     "./Clemens/example_sym_2/ENSG00000138316_ADAMTS14.fasta",
#     "./Clemens/toy_example_GTR+G4/toy_example_ntaxa_7_run_1-alignment.phy",
#     "./test/toy_example_only_msa/toy_example_ntaxa_7_run_5-alignment.phy",
# ]

if __name__ == "__main__":
    for test_case in test_cases_msa_tree_model:
        print(f"Processing {test_case}")
        process_with_msa_tree_model(test_case)
    # for test_case in test_case_without_model:
    #     print(f"Processing {test_case}")
    #     process_without_model(test_case)
    # for test_case in test_cases_msa:
    #     print(f"Processing {test_case}")
    #     process_without_model(test_case)
