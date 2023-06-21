import numpy as np
import regex as re
import pandas as pd
from ete3 import Tree
import scipy
import scipy.linalg
import scipy.stats as st
import subprocess
from pathlib import Path
import os
from rich import print

from satute_repository import (
    parse_rate_and_frequencies_and_create_model_files,
    parse_rate_matrices,
    parse_output_state_frequencies,
)
from satute_test_statistic import calculate_test_statistic
from satute_exception import InputArgumentsError

def remove_filename(path):
    parts = path.split("/")
    nameFILE = parts[-2]
    parts.pop()
    pathFOLDER = "/".join(parts) + "/"  # Path where the we create folders
    return pathFOLDER

"""## INTERNAL NODES AND LEAVES"""


def node_type(T):
    leaves = []
    for i in T.get_leaves():
        leaves.append(i.name)

    internal_nodes = []
    for node in T.traverse("levelorder"):
        if (node not in T.get_leaves()) and (node in T.get_tree_root()):
            internal_nodes.append(node.name)

    return leaves, internal_nodes


"""## MODIFYING SEQUENCE FILE

#### We modified the nodes names in order to have indistinguishable names. We need to modify as well the sequence file.
"""


def modify_seqfile(path, leaves, option):
    filesequence = path

    if option == 1:  # assumption: msa in  phylip format
        with open(filesequence, "r+") as f:
            with open(path + ".txt", "w") as writer:
                lines = f.readlines()
                writer.write(lines[0])
                for l in lines:
                    for i in leaves:
                        if l[0 : len(i) - 1] == i[0 : len(i) - 1]:
                            lmod = l[: len(i) - 1] + "*" + l[len(i) - 1 :]
                            writer.write(lmod)

    elif option == 2:  # assumption:  msa in fasta file
        with open(filesequence, "r+") as f:
            with open(path + ".txt", "w") as writer:
                lines = f.readlines()
                for line in lines:
                    if line[0] == ">":
                        for i in leaves:
                            if line[1:-1] == i[0 : len(i) - 1]:
                                writer.write(">" + i + "\n")
                    else:
                        writer.write(line)

    else:
        print("FILE FORMAT OPTION NOT VALID")


"""## COMPUTING BRANCHES LENGTHS """


def branch_lengths(T):
    vector_branches = []
    vector_distances = []
    internal_nodes = []

    for node in T.traverse("levelorder"):
        if (node not in T.get_leaves()) and (node in T.get_tree_root()):
            internal_nodes.append(node.name)

    for node in T.traverse("levelorder"):  # First internal branches.
        children = node.get_children()

        for child in children:
            if child.name in internal_nodes:
                vector_branches.append(node.name + "-" + child.name)
                vector_distances.append(T.get_distance(node, child))

    for node in T.traverse("levelorder"):  # Same for external branches.
        children = node.get_children()

        for child in children:
            if child.name not in internal_nodes:
                vector_branches.append(node.name + "-" + child.name)
                vector_distances.append(T.get_distance(node, child))

    return vector_branches, vector_distances


"""## NUMBER OF NODES AND SITES"""


def number_nodes_sites(path):
    full_posterior_probabilities_file = f"{path}.state"
    tree_posterior_probabilities_table = parse_output_state_frequencies(
        full_posterior_probabilities_file
    )
    nodes_number = len(tree_posterior_probabilities_table["Node"].unique())
    nucleotides_sites = len(tree_posterior_probabilities_table["Site"].unique())
    return nodes_number, nucleotides_sites


"""## SEPARATING CLADES"""


def clades(T, t, newick_format, internal_nodes, leaves):
    root = T.get_tree_root()

    parent_nodes = []
    for i in internal_nodes:
        node = T.search_nodes(name=i)[0]
        parent = node.up
        parent_nodes.append(parent.name)

    clades1 = []
    clades2 = []

    for i in internal_nodes:
        t1 = t
        cont = 0
        length = len(i)

        for j in range(len(t1)):
            if t1[j : j + length] == i:
                for k in reversed(range(len(t1[0:j]))):
                    if t1[k] == ")":
                        cont += 1
                    elif t1[k] == "(":
                        cont -= 1
                        aux = k

                    if cont == 0:
                        clade1 = t1[k : j + length] + ";"
                        clades1.append(clade1)
                        break

        root_clade2 = parent_nodes[
            internal_nodes.index(i)
        ]  # parent node of current node that will be the root of clade 2

        if root_clade2 == root.name:
            clade1 = clade1.replace(";", "")
            clade2 = t1.replace(clade1, "")
            for j in range(len(clade2)):
                if clade2[j : j + 2] == ",:":
                    split = clade2.split(",:", 1)[1]
                    if split.find(",(") > 0:
                        clade2 = clade2.replace(clade2[j : j + split.find(",(")], "", 1)
                    else:
                        clade2 = clade2.replace(clade2[j : j + split.find(")")], "", 1)

        else:
            r = root_clade2

            # FIRST WE NEED TO SAVE THE NODES FROM ROOT TO END AND THEIR CORRESPONDING CHILDREN

            auxiliar = []
            auxiliar_parent = []
            vect_fromroottoend = []

            r1 = r
            current = i
            while r1 != root.name:
                vect_fromroottoend.append(r)
                children = T.search_nodes(name=r)[0].get_children()

                for k in range(len(children)):
                    if children[k].name != current:
                        if children[k].name in leaves:
                            split = t1.split(children[k].name + ":", 1)[1]
                            if split[0 : split.find(")")] < split[0 : split.find(",")]:
                                auxiliar.append(
                                    children[k].name + ":" + split[0 : split.find(")")]
                                )
                            else:
                                auxiliar.append(
                                    children[k].name + ":" + split[0 : split.find(",")]
                                )
                            auxiliar_parent.append(r)
                        else:
                            node = T.search_nodes(name=children[k].name)[0]
                            node = node.write(format=newick_format)
                            node = node.replace(";", "")
                            auxiliar.append(node)
                            auxiliar_parent.append(r)

                r1 = r
                if r1 != root.name:
                    r = T.search_nodes(name=r)[0].up.name
                    current = r1

            # NOW WE SAVE DISTANCES

            distances_fromroottoend = [""]
            for j in range(1, len(vect_fromroottoend)):
                split = t1.split(vect_fromroottoend[j - 1] + ":", 1)[1]
                if split[0 : split.find(")")] < split[0 : split.find(",")]:
                    distances_fromroottoend.append(split[0 : split.find(")")])
                else:
                    distances_fromroottoend.append(split[0 : split.find(",")])

            # LET'S CONSTRUCT THE NEW CLADE

            for j in reversed(vect_fromroottoend):
                index = vect_fromroottoend.index(j)
                if j != root_clade2:
                    if auxiliar_parent.count(j) == 2:
                        attach = []
                        for k in range(len(auxiliar_parent)):
                            if auxiliar_parent[k] == j:
                                attach.append(auxiliar[k])
                        clade2 = (
                            "("
                            + attach[0]
                            + ","
                            + attach[1]
                            + ")"
                            + j
                            + ":"
                            + distances_fromroottoend[index]
                        )
                    else:
                        clade2 = (
                            "("
                            + auxiliar[index]
                            + ","
                            + clade2
                            + ")"
                            + j
                            + ":"
                            + distances_fromroottoend[index]
                        )
                else:
                    clade2 = (
                        "(" + auxiliar[index] + "," + clade2 + ")" + root_clade2 + ";"
                    )

        clades2.append(clade2)

    parent_leaves = []
    for i in leaves:
        node = T.search_nodes(name=i)[0]
        parent = node.up
        parent_leaves.append(parent.name)

    for i in range(len(leaves)):
        clade1 = (
            "(copy_" + leaves[i] + ":0.0001000000," + leaves[i] + ":0.0001000000)ROOT;"
        )
        clades1.append(clade1)

        root_clade2 = parent_leaves[
            i
        ]  # parent node of leaf that will be the root of clade 2

        if root_clade2 == root.name:
            split = t.split(leaves[i], 1)[1]
            clade2 = t.replace(leaves[i], "")
            clade2 = clade2.replace(split[0 : split.find(",") + 1], "", 1)

        else:
            r = root_clade2

            # FIRST WE NEED TO SAVE THE NODES FROM ROOT TO END AND THEIR CORRESPONDING CHILDREN

            auxiliar = []
            auxiliar_parent = []
            vect_fromroottoend = []

            r1 = r
            current = leaves[i]
            while r1 != root.name:
                vect_fromroottoend.append(r)
                children = T.search_nodes(name=r)[0].get_children()

                for k in range(len(children)):
                    if children[k].name != current:
                        if children[k].name in leaves:
                            split = t1.split(children[k].name + ":", 1)[1]
                            if split[0 : split.find(")")] < split[0 : split.find(",")]:
                                auxiliar.append(
                                    children[k].name + ":" + split[0 : split.find(")")]
                                )
                            else:
                                auxiliar.append(
                                    children[k].name + ":" + split[0 : split.find(",")]
                                )
                            auxiliar_parent.append(r)
                        else:
                            node = T.search_nodes(name=children[k].name)[0]
                            node = node.write(format=newick_format)
                            node = node.replace(";", "")
                            auxiliar.append(node)
                            auxiliar_parent.append(r)

                r1 = r
                if r1 != root.name:
                    r = T.search_nodes(name=r)[0].up.name
                    current = r1

            # NOW WE SAVE DISTANCES

            distances_fromroottoend = [""]
            for j in range(1, len(vect_fromroottoend)):
                split = t1.split(vect_fromroottoend[j - 1] + ":", 1)[1]
                if split[0 : split.find(")")] < split[0 : split.find(",")]:
                    distances_fromroottoend.append(split[0 : split.find(")")])
                else:
                    distances_fromroottoend.append(split[0 : split.find(",")])

            # LET'S CONSTRUCT THE NEW CLADE

            for j in reversed(vect_fromroottoend):
                index = vect_fromroottoend.index(j)
                if j != root_clade2:
                    if auxiliar_parent.count(j) == 2:
                        attach = []
                        for k in range(len(auxiliar_parent)):
                            if auxiliar_parent[k] == j:
                                attach.append(auxiliar[k])
                        clade2 = (
                            "("
                            + attach[0]
                            + ","
                            + attach[1]
                            + ")"
                            + j
                            + ":"
                            + distances_fromroottoend[index]
                        )
                    else:
                        clade2 = (
                            "("
                            + auxiliar[index]
                            + ","
                            + clade2
                            + ")"
                            + j
                            + ":"
                            + distances_fromroottoend[index]
                        )
                else:
                    clade2 = (
                        "(" + auxiliar[index] + "," + clade2 + ")" + root_clade2 + ";"
                    )
        clades2.append(clade2)

    return clades1, clades2


def run_iqtree_for_each_clade(path_folder, number_rates, chosen_rate, iqtree_path):
    """Prepares necessary information and runs the IQ-TREE."""
    model_and_frequency = ""
    path_new_folder = ""

    # Condition to define path and model
    if number_rates > 1:
        path_new_folder = os.path.join(
            path_folder, "subsequences", f"subseq{chosen_rate}", "clades"
        )
        model_and_frequency_file = os.path.join(
            path_folder, "subsequences", f"subseq{chosen_rate}", "model.txt"
        )
    else:
        path_new_folder = os.path.join(path_folder, "clades")
        model_and_frequency_file = os.path.join(path_folder, "model.txt")

    # Check if model and frequency file exists
    if not os.path.isfile(model_and_frequency_file):
        raise FileNotFoundError(
            f"The model and frequency file '{model_and_frequency_file}' does not exist."
        )

    # Read model and frequency
    with open(model_and_frequency_file, "r") as toModel:
        model_and_frequency = toModel.readline().strip()

    # Check if path_new_folder exists and is a directory
    if not os.path.isdir(path_new_folder):
        raise NotADirectoryError(
            f"The path '{path_new_folder}' does not exist or is not a directory."
        )

    # Iterate over each clade directory
    for clade_dir in Path(path_new_folder).iterdir():
        if clade_dir.is_dir():
            # Command to be executed for each clade
            cmd = [
                iqtree_path,
                "-s",
                "sequence.txt",
                "-te",
                "tree.txt",
                "-m",
                model_and_frequency,
                "-asr",
                "-blfix",
                "-o",
                "FOO",
                "-pre",
                "output",
                "-redo",
                "-quiet",
            ]

            # Check if sequence and tree files exist in the clade directory
            if not os.path.isfile(
                os.path.join(clade_dir, "sequence.txt")
            ) or not os.path.isfile(os.path.join(clade_dir, "tree.txt")):
                raise FileNotFoundError(
                    "Either sequence.txt or tree.txt file does not exist in directory: "
                    + str(clade_dir)
                )

            # Run the command
            result = subprocess.run(cmd, cwd=clade_dir)

            # Check if the command was successful
            if result.returncode != 0:
                raise RuntimeError(
                    f"The command '{' '.join(cmd)}' failed with return code: {result.returncode}"
                )


import os
from multiprocessing import Pool
from functools import partial

def execute_iqtree(sub_dir, iqtree_path, model_and_frequency):
    # Command to be executed for each clade
    cmd = [
        iqtree_path,
        "-s",
        "sequence.txt",
        "-te",
        "tree.txt",
        "-m",
        model_and_frequency,
        "-asr",
        "-blfix",
        "-o",
        "FOO",
        "-pre",
        "output",
        "-redo",
        "-quiet",
    ]

    # Check if sequence and tree files exist in the subdirectory
    if not os.path.isfile(os.path.join(sub_dir, "sequence.txt")) or not os.path.isfile(
        os.path.join(sub_dir, "tree.txt")
    ):
        raise FileNotFoundError(
            f"Either sequence.txt or tree.txt file does not exist in directory: {sub_dir}"
        )

    # Run the IQ-TREE command
    result = subprocess.run(cmd, cwd=sub_dir)

    # Check if the command was successful
    if result.returncode != 0:
        raise RuntimeError(
            f"The command '{' '.join(cmd)}' failed with return code: {result.returncode}"
        )

def run_iqtree_for_each_clade_parallel(path_folder, number_rates, chosen_rate, iqtree_path):
    """
    Run IQ-TREE for each clade directory in parallel.

    Args:
        path_folder (str): Root folder path.
        number_rates (int): Number of rates.
        chosen_rate (int): Chosen rate.
        iqtree_path (str): Path to the IQ-TREE executable.

    Raises:
        FileNotFoundError: If model and frequency file does not exist.
        NotADirectoryError: If the clade directory does not exist or is not a directory.
        FileNotFoundError: If sequence.txt or tree.txt file does not exist in the clade directory.
        RuntimeError: If IQ-TREE command fails.

    """
    model_and_frequency = ""
    clade_dir = ""

    # Determine path and model based on number of rates
    if number_rates > 1:
        clade_dir = os.path.join(
            path_folder, "subsequences", f"subseq{chosen_rate}", "clades"
        )
        model_and_frequency_file = os.path.join(
            path_folder, "subsequences", f"subseq{chosen_rate}", "model.txt"
        )
    else:
        clade_dir = os.path.join(path_folder, "clades")
        model_and_frequency_file = os.path.join(path_folder, "model.txt")

    # Check if model and frequency file exists
    if not os.path.isfile(model_and_frequency_file):
        raise FileNotFoundError(
            f"The model and frequency file '{model_and_frequency_file}' does not exist."
        )

    # Read model and frequency
    with open(model_and_frequency_file, "r") as toModel:
        model_and_frequency = toModel.readline().strip()

    # Check if clade directory exists and is a directory
    if not os.path.isdir(clade_dir):
        raise NotADirectoryError(
            f"The clade directory '{clade_dir}' does not exist or is not a directory."
        )

    # Create a list of clade directories
    clade_dirs = [sub_dir for sub_dir in Path(clade_dir).iterdir() if sub_dir.is_dir()]

    # Create a partial function with fixed arguments for execute_iqtree
    execute_iqtree_partial = partial(execute_iqtree, iqtree_path=iqtree_path, model_and_frequency=model_and_frequency)

    # Create a multiprocessing pool and map the execute_iqtree_partial function to clade_dirs
    with Pool() as pool:
        pool.map(execute_iqtree_partial, clade_dirs)



"""## DIVIDING SEQUENCE FILE INTO SUBFILES DEPENDING ON WHICH RATE IS MOST PROBABLE"""


def subsequences(T, path, epsilon, number_rates, option):
    filesequence = path + ".txt"

    fileprob = path + ".siteprob"

    pathFolder = remove_filename(path)
    with open(fileprob, "r+") as f:
        with open(pathFolder + "prob.csv", "w") as writer:
            lines = f.readlines()
            for i in range(len(lines)):
                writer.write(lines[i])

    df = pd.read_csv(pathFolder + "prob.csv", sep="\t", engine="python")
    site_rate = []
    length = len(df.index)
    for i in range(length):
        probs = []
        for j in range(number_rates):
            probs.append(df.iloc[i][j + 1])
        max_value = max(probs)
        if max_value < 1 / number_rates + epsilon:
            site_rate.append(0)
        else:
            site_rate.append(probs.index(max_value) + 1)
    numbersitesperrate = []
    for j in range(number_rates):
        numbersitesperrate.append(site_rate.count(j + 1))

    if option == 1:
        for i in range(number_rates):
            os.makedirs(pathFolder + "subsequences/subseq" + str(i + 1), exist_ok=True)
            fseq = open(
                pathFolder + "subsequences/subseq" + str(i + 1) + "/sequence.txt", "w+"
            )
            with open(filesequence, "r+") as f:
                lines = f.readlines()
                l = lines[0]
                fseq.write(
                    l[0 : l.find(" ")] + " " + str(site_rate.count(i + 1)) + "\n"
                )
                for j in range(1, len(lines)):
                    line = lines[j]
                    fseq.write(line.split(" ", 1)[0])
                    seq = line.split(" ", 1)[1]
                    for k in range(len(seq)):
                        if seq[k] != " ":
                            index = k
                            break
                    seq = seq[index:-1]
                    # cassius
                    fseq.write((k + 1) * " ")
                    for k in range(len(site_rate)):
                        if site_rate[k] == i + 1:
                            fseq.write(seq[k])
                    fseq.write("\n")
    elif option == 2:
        leaves = []
        for i in T.get_leaves():
            leaves.append(i.name)
        for i in range(number_rates + 1):
            os.makedirs(pathFolder + "subsequences/subseq" + str(i + 1), exist_ok=True)
            fseq = open(
                pathFolder + "subsequences/subseq" + str(i + 1) + "/sequence.txt", "w+"
            )
            with open(filesequence, "r+") as f:
                lines = f.readlines()
                seq = ""
                for j in range(len(lines)):
                    if lines[j][1:-1] in leaves:
                        if len(seq) > 0:
                            for k in range(len(site_rate)):
                                if site_rate[k] == i + 1:
                                    fseq.write(seq[k])
                            fseq.write("\n")
                        fseq.write(lines[j])
                        seq = ""
                    else:
                        seq = seq + lines[j][0 : lines[j].find(" ")]
            for k in range(len(site_rate)):
                if site_rate[k] == i + 1:
                    fseq.write(seq[k])
    else:
        print("FILE FORMAT OPTION NOT VALID")

    return numbersitesperrate


"""## SAVING MUTATION RATES (in case of the Gamma model)"""


def save_rates(path, number_rates):
    rates = []
    filne = path + ".iqtree"
    with open(filne, "r+") as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if " Category  Relative_rate  Proportion" in line:
                if lines[i + 1][2:3] == 0:
                    for j in range(number_rates):
                        rates.append(float(lines[i + j + 2][12:20]))
                else:
                    for j in range(number_rates):
                        rates.append(float(lines[i + j + 1][12:20]))
    return rates


"""## MODIFY BRANCHES LENGTHS (if it corresponds), ADD FOO AND SAVE CLADES"""


def save_clades(path, number_rates, clades1, clades2, newick_format, rates):
    pathFolder = remove_filename(path)
    if number_rates == 1:
        for i in range(len(clades1)):
            os.makedirs(
                pathFolder + "clades/Branch" + str(i) + "_clade1/", exist_ok=True
            )
            os.makedirs(
                pathFolder + "clades/Branch" + str(i) + "_clade2/", exist_ok=True
            )

            #f"(FOO:1((A:1,B:1):1,C:1));"
            # f"FOO:0.00000000010,({newickformat.remove(";")});"

            clades1[i] = (
                "(FOO:0.00000000010,"
                + clades1[i][0 : clades1[i].rfind(")")]
                + ")"
                + clades1[i][clades1[i].rfind(")") : -1]
                + ";"
            )
            clades2[i] = (
                "(FOO:0.00000000010,"
                + clades2[i][0 : clades2[i].rfind(")")]
                + ")"
                + clades2[i][clades2[i].rfind(")") : -1]
                + ";"
            )

            f1 = open(pathFolder + "clades/Branch" + str(i) + "_clade1/tree.txt", "w")
            f1.write(clades1[i])

            f2 = open(pathFolder + "clades/Branch" + str(i) + "_clade2/tree.txt", "w")
            f2.write(clades2[i])

    else:
        for j in range(number_rates):
            for i in range(len(clades1)):
                cl1 = clades1[i]
                cl2 = clades2[i]

                os.makedirs(
                    pathFolder
                    + "subsequences/subseq"
                    + str(j + 1)
                    + "/clades/Branch"
                    + str(i)
                    + "_clade1/",
                    exist_ok=True,
                )
                os.makedirs(
                    pathFolder
                    + "subsequences/subseq"
                    + str(j + 1)
                    + "/clades/Branch"
                    + str(i)
                    + "_clade2/",
                    exist_ok=True,
                )

                C1 = Tree(cl1, format=newick_format)
                r = C1.get_tree_root()
                for node in C1.traverse("levelorder"):
                    node.dist = rates[j] * node.dist
                c = C1.write(format=newick_format)
                cl1 = c[0:-1] + r.name + ";"

                C2 = Tree(cl2, format=newick_format)
                r = C2.get_tree_root()
                for node in C2.traverse("levelorder"):
                    node.dist = rates[j] * node.dist
                c = C2.write(format=newick_format)
                cl2 = c[0:-1] + r.name + ";"

                cl1 = (
                    "(FOO:0.00000000010,"
                    + cl1[0 : cl1.rfind(")")]
                    + ")"
                    + cl1[cl1.rfind(")") : -1]
                    + ";"
                )
                cl2 = (
                    "(FOO:0.00000000010,"
                    + cl2[0 : cl2.rfind(")")]
                    + ")"
                    + cl2[cl2.rfind(")") : -1]
                    + ";"
                )

                f1 = open(
                    pathFolder
                    + "subsequences/subseq"
                    + str(j + 1)
                    + "/clades/Branch"
                    + str(i)
                    + "_clade1/tree.txt",
                    "w",
                )
                f1.write(cl1)

                f2 = open(
                    pathFolder
                    + "subsequences/subseq"
                    + str(j + 1)
                    + "/clades/Branch"
                    + str(i)
                    + "_clade2/tree.txt",
                    "w",
                )
                f2.write(cl2)


"""## MODIFYING BRANCHES LENGTHS IN FULL TREES"""


def modify_fulltree(path, T, rates, newick_format):
    for j in range(rates):
        r = T.get_tree_root()
        for node in T.traverse("levelorder"):
            node.dist = rates[j] * node.dist
        t = T.write(format=newick_format)
        t = t[0:-1] + r.name + ";"

        f = open(path + "subsequences/subseq" + str(j + 1) + "/tree.txt", "w")
        f.write(t)


"""## CREATING FOR EACH CLADE THE NUCLEOTIDE SEQUENCE FILE"""


def sequences_clades(
    path,
    number_rates,
    nodes_number,
    nucleotides_sites,
    clades1,
    clades2,
    option,
    newick_format,
    internal_nodes,
    numbersitesperrate,
):
    
    
    pathFolder = remove_filename(path)

    if number_rates == 1:
        if option == 1:
            filesequence = path + ".txt"

            for i in range(len(clades1)):
                C1 = Tree(clades1[i], format=newick_format)

                f1 = open(
                    pathFolder + "clades/Branch" + str(i) + "_clade1/sequence.txt", "w"
                )

                leaves = []
                for k in C1.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, "r+") as f:
                    lines = f.readlines()
                    f1.write(str(len(leaves)) + lines[0][lines[0].find(" ") : -1])
                    f1.write("\n")
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[0 : line.find(" ")] in leaves:
                            f1.write(line)

                if i >= len(
                    internal_nodes
                ):  # if clade is only a leaf, we need to rewrite it in the sequence file
                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0 : line.find(" ")] in leaves:
                                f1.write("copy_" + line)

                f1.write(
                    "FOO  "
                    + "N" * (int(lines[0][lines[0].find(" ") + 1 : -1]) - 1)
                    + "A"
                )

            for i in range(len(clades2)):
                C2 = Tree(clades2[i], format=newick_format)

                f2 = open(
                    pathFolder + "clades/Branch" + str(i) + "_clade2/sequence.txt", "w"
                )

                leaves = []
                for k in C2.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, "r+") as f:
                    lines = f.readlines()
                    f2.write(str(len(leaves)) + lines[0][lines[0].find(" ") : -1])
                    f2.write("\n")
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[0 : line.find(" ")] in leaves:
                            f2.write(line)

                f2.write(
                    "FOO  "
                    + "N" * (int(lines[0][lines[0].find(" ") + 1 : -1]) - 1)
                    + "A"
                )

        elif option == 2:
            filesequence = path + ".txt"

            for i in range(len(clades1)):
                C1 = Tree(clades1[i], format=newick_format)

                f1 = open(
                    pathFolder + "clades/Branch" + str(i) + "_clade1/sequence.txt", "w"
                )

                leaves = []
                for k in C1.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, "r+") as f:
                    lines = f.readlines()
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[1:-1] in leaves:
                            f1.write(line)
                            for k in range(j + 1, len(lines)):
                                line = lines[k]
                                if line[0] == ">":
                                    break
                                else:
                                    f1.write(line)

                if i >= len(
                    internal_nodes
                ):  # if clade is only a leaf, we need to rewrite it in the sequence file
                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f1.write(line[0] + "copy_" + line[1:-1] + "\n")
                                for k in range(j + 1, len(lines)):
                                    line = lines[k]
                                    if line[0] == ">":
                                        break
                                    else:
                                        f1.write(line)

                f1.write(">FOO\n" + "N" * (nucleotides_sites - 1) + "A")

            for i in range(len(clades2)):
                C2 = Tree(clades2[i], format=newick_format)

                f2 = open(
                    pathFolder + "clades/Branch" + str(i) + "_clade2/sequence.txt", "w"
                )

                leaves = []
                for k in C2.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, "r+") as f:
                    lines = f.readlines()
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[1:-1] in leaves:
                            f2.write(line)
                            for k in range(j + 1, len(lines)):
                                line = lines[k]
                                if line[0] == ">":
                                    break
                                else:
                                    f2.write(line)

                f2.write(">FOO\n" + "N" * (nucleotides_sites - 1) + "A")

        else:
            print("FILE FORMAT OPTION NOT VALID")

    else:
        if option == 1:
            for r in range(number_rates):
                filesequence = (
                    pathFolder + "subsequences/subseq" + str(r + 1) + "/sequence.txt"
                )

                for i in range(len(clades1)):
                    C1 = Tree(clades1[i], format=newick_format)

                    f1 = open(
                        pathFolder
                        + "subsequences/subseq"
                        + str(r + 1)
                        + "/clades/Branch"
                        + str(i)
                        + "_clade1/sequence.txt",
                        "w",
                    )

                    leaves = []
                    for k in C1.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        f1.write(
                            str(len(leaves) + 1) + lines[0][lines[0].find(" ") : -1]
                        )
                        f1.write("\n")
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0 : line.find(" ")] in leaves:
                                f1.write(line)

                    if i >= len(
                        internal_nodes
                    ):  # if clade is only a leaf, we need to rewrite it in the sequence file
                        with open(filesequence, "r+") as f:
                            lines = f.readlines()
                            for j in range(0, len(lines)):
                                line = lines[j]
                                if line[0 : line.find(" ")] in leaves:
                                    f1.write("copy_" + line)

                    f1.write("FOO  " + "N" * (numbersitesperrate[r] - 1) + "A")

                for i in range(len(clades2)):
                    C2 = Tree(clades2[i], format=newick_format)

                    f2 = open(
                        pathFolder
                        + "subsequences/subseq"
                        + str(r + 1)
                        + "/clades/Branch"
                        + str(i)
                        + "_clade2/sequence.txt",
                        "w",
                    )

                    leaves = []
                    for k in C2.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        f2.write(
                            str(len(leaves) + 1) + lines[0][lines[0].find(" ") : -1]
                        )
                        f2.write("\n")
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0 : line.find(" ")] in leaves:
                                f2.write(line)

                    f2.write("FOO  " + "N" * (numbersitesperrate[r] - 1) + "A")

        elif option == 2:
            for r in range(number_rates):
                filesequence = (
                    pathFolder + "subsequences/subseq" + str(r + 1) + "/sequence.txt"
                )

                for i in range(len(clades1)):
                    C1 = Tree(clades1[i], format=newick_format)

                    f1 = open(
                        pathFolder
                        + "subsequences/subseq"
                        + str(r + 1)
                        + "/clades/Branch"
                        + str(i)
                        + "_clade1/sequence.txt",
                        "w",
                    )

                    leaves = []
                    for k in C1.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f1.write(line)
                                for k in range(j + 1, len(lines)):
                                    line = lines[k]
                                    if line[0] == ">":
                                        break
                                    else:
                                        f1.write(line)
                                f1.write("\n")

                    if i >= len(
                        internal_nodes
                    ):  # if clade is only a leaf, we need to rewrite it in the sequence file
                        with open(filesequence, "r+") as f:
                            lines = f.readlines()
                            for j in range(0, len(lines)):
                                line = lines[j]
                                if line[1:-1] in leaves:
                                    f1.write(line[0] + "copy_" + line[1:-1] + "\n")
                                    for k in range(j + 1, len(lines)):
                                        line = lines[k]
                                        if line[0] == ">":
                                            break
                                        else:
                                            f1.write(line)
                                    f1.write("\n")

                    f1.write(">FOO\n" + "N" * (numbersitesperrate[r] - 1) + "A")

                for i in range(len(clades2)):
                    C2 = Tree(clades2[i], format=newick_format)

                    f2 = open(
                        pathFolder
                        + "subsequences/subseq"
                        + str(r + 1)
                        + "/clades/Branch"
                        + str(i)
                        + "_clade2/sequence.txt",
                        "w",
                    )

                    leaves = []
                    for k in C2.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, "r+") as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f2.write(line)
                                for k in range(j + 1, len(lines)):
                                    line = lines[k]
                                    if line[0] == ">":
                                        break
                                    else:
                                        f2.write(line)
                                f2.write("\n")

                    f2.write(">FOO\n" + "N" * (numbersitesperrate[r] - 1) + "A")

        else:
            print("FILE FORMAT OPTION NOT VALID")


def compare_arrays(array1, array2):
    """Compare two numpy arrays for equality and print a message."""
    if np.array_equal(array1, array2):
        print("The arrays are identical! 😎")
    else:
        print("The arrays are not identical. 😞")


"""## SPECTRAL DECOMPOSITION OF THE RATE MATRIX"""

def spectral_decomposition(n, path):
    rate_matrix, psi_matrix = parse_rate_matrices(n, path)

    """ 
        Then psi_matrix := Diag(pi). Recall that matrix Q is reversible iff M:= psi_matrix^1/2 x Q x psi_matrix^{-1/2} is symmetric.
        For a real symmetric matrix M, its eigenvectors can be chosen to be an orthonormal basis of R^n 
    """
    M = scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2) @ rate_matrix
    M = M @ scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2)

    """ eigendecomposition of matrix M"""
    lamb, w = np.linalg.eig(M)  # Compute the eigenvalues and eigenvectors.
    idx = lamb.argsort()[::-1]  # Order from large to small.
    lamb = lamb[idx] # Order the eigenvalues (large to small).
    w = w[:, idx] # Order the eigenvectors according to the eigenvalues"""

    # the first one should be the eigenvalue 0 in lamb, why are we doing the following?
    lamb_nozero = []  # list of eigenvalues without 0
    for i in lamb:
        if i > 0.00999 or i < -0.00999:
            lamb_nozero.append(i)
        
    max_lambda = max(lamb_nozero)  # dominant non-zero eigenvalue 
    # get the indices of the dominant non-zero eigenvalue in lamb taking numerical inaccuracies into account and identical values
    index = []
    for i in range(len(lamb)):
        lambda_it = lamb[i]
        if abs(lambda_it - max_lambda) < 0.01:
            index.append(i)

    multiplicity = len(index) # multiplicity of the dominant non-zero eigenvalue
    array_eigenvectors = [] # list of right eigenvectors for the dominant non-zero eigenvalue
    for i in range(multiplicity):
        # calculate the right eigenvectors for the dominant non-zero eigenvalue
        v1 = (
            scipy.linalg.fractional_matrix_power(psi_matrix, -1 / 2)
            @ w[:, index[i]]
            )
        array_eigenvectors.append(v1)
        """# calculate the left eigenvectors for the dominant non-zero eigenvalue
        h1 = (
            scipy.linalg.fractional_matrix_power(psi_matrix, +1 / 2)
            @ w[:, index[i]]
            ) """
        
    return array_eigenvectors, multiplicity


""" Read newick string from convert to ete3 satute format """


def convert_newick_to_satute_ete3_format(t, newick_format):
    #print(newick_format)
    T = Tree(t, format=newick_format)
    print(t)
    print(T.write(format=newick_format))

    for node in T.traverse("levelorder"):
        l = len(node.name)
        for i in range(len(t)):
            if t[i : i + l] == str(node.name) and (t[i + l] == ";" or t[i + l] == ":"):
                t = t[: i + l] + "*" + t[i + l :]

    T = Tree(t, format=newick_format)

    root = T.get_tree_root()
    root_children = root.get_children()
    leaves = T.get_leaves()

    root_children_leaves = []
    for i in root_children:
        if i in leaves:
            root_children_leaves.append(i)

    if len(root_children_leaves) >= 2:
        for i in root_children_leaves:
            for j in range(len(t)):
                if t[j : j + len(i.name)] == i.name:
                    cont = 0
                    split0 = t.split(i.name, 1)[0]
                    split1 = t.split(i.name, 1)[1]
                    for k in range(len(split0)):
                        if split0[k] == "(":
                            cont += 1
                    if cont > 1:
                        t = t.replace(i.name, "")
                        if len(split1[0 : split1.find(")")]) < len(
                            split1[0 : split1.find(",")]
                        ):
                            t = t.replace(split1[0 : split1.find(")")], "")
                            t = (
                                t[0]
                                + str(i.name)
                                + split1[0 : split1.find(")")]
                                + ","
                                + t[1:-1]
                            )
                        else:
                            t = t.replace(split1[0 : split1.find(",")], "")
                            t = (
                                t[0]
                                + str(i.name)
                                + split1[0 : split1.find(",")]
                                + ","
                                + t[1:-1]
                            )

    for i in range(len(t)):
        if t[i : i + 2] == ",)":
            t = t.replace(t[i : i + 2], ")")

    T = Tree(t, format=newick_format)

    return t, T


def guess_msa_file_format(file_path):
    """
    Guess the format of a file based on its first line.

    :param file_path: Path to the file
    :return: The guessed file format or None if the format could not be guessed
    """
    with open(file_path, "r") as file:
        first_line = file.readline().strip()

    # PHYLIP files typically start with two integers (number of species and number of characters)
    line_parts = first_line.split()
    if len(line_parts) == 2 and line_parts[0].isdigit() and line_parts[1].isdigit():
        return 1  ##'PHYLIP'
    # FASTA files typically start with a '>' character
    elif first_line.startswith(">"):
        return 2  #'FASTA'
    #else:
    #    raise InputArgumentsError("Wrong msa file format!")



def write_results_and_newick_tree(
    results_list, newick_string, path_folder, chosen_rate, c_sTwoSequence, T
):
    """
    This function writes the saturation branches data to a file and then appends the saturation information
    newick string to the same file.

    :param results_list: List of results data to be written to file.
    :param newick_string: Newick formatted string representing the tree.
    :param path_folder: Path to the folder where results will be saved.
    :param chosen_rate: Chosen rate parameter to be included in the filename.
    :param c_sTwoSequence: Saturation coherence between two sequences.
    :param T: ETE Tree instance.
    """

    # Convert the results_list into a pandas dataframe
    saturation_branches_dataframe = pd.DataFrame(results_list)

    # Save the dataframe as a tab-separated CSV file
    saturation_branches_dataframe.to_csv(
        f"{path_folder}/resultsRate{chosen_rate}.satute.csv",
        header=True,
        index=None,
        sep="\t",
        mode="a",
    )



    # Generate a newick string with saturation information
    saturation_information_newick_string = map_values_to_newick(
        results_list, newick_string
    )

    # Open the results file in append mode
    with open(
        f"{path_folder}/resultsRate{chosen_rate}.satute.csv", "a"
    ) as satute_result_file:
        # Write additional information to the file
        satute_result_file.write(
            "\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is  {:6.4f}".format(
                c_sTwoSequence
            )
        )

        satute_result_file.write(
            "\n\nFor better reference, this is the reconstructed tree topology :\n\n"
        )

        # Write the ascii representation of the tree to the file
        satute_result_file.write(
            T.copy("newick").get_ascii(attributes=["name", "label", "distance"])
        )

        # Tree without saturation values            
        satute_result_file.write(f"\n\n Tree with saturation values: {newick_string}")
        # Write the saturation information newick string to the file
        satute_result_file.write(f"\n\n Tree with saturation values: {saturation_information_newick_string}")


""" # MAIN FUNCTION FOR THE SATURATION TEST """
def saturation_test_cli(
    pathDATA,
    newick_string,
    pathIQTREE,
    dimension=4,
    number_rates=4,
    chosen_rate=str(4),
    alpha=0.01,
    newick_format=1,
    epsilon=0.01,
    model="GTR",
):
    """
    :param pathDATA: str
        A string representing the path to the data file. This file should be in the appropriate format for the analysis to be performed.

    :param pathIQTREE: str
        A string representing the path to the IQ-TREE executable file.

    :param dimension: int, default = 4
        The dimension of the data. This refers to the number of distinct symbols in the sequence data (e.g., for DNA, dimension would typically be 4 corresponding to A, C, G, T).

    :param number_rates: int, default = 4
        The number of rate categories to use in the analysis.

    :param chosen_rate: str, default = '4'
        The rate category to be used for the analysis. The rates are numbered starting from 1.

    :param alpha: float, default = 0.01
        The significance level of the test.

    :param newick_format: int, default = 1
        The format of the output tree file. The default format (1) is the standard Newick format.

    :param epsilon: float, default = 0.01
        A small positive number used as a tolerance in numerical calculations.

    """
    results_list = []

    """Check type if the alignment is Phylip or fasta format"""
    option = guess_msa_file_format(pathDATA)

    """ PATH  to input directory Need to be changed """
    pathFOLDER = remove_filename(str(pathDATA))

    """ data structure for phylogenetic tree"""
    t, T = convert_newick_to_satute_ete3_format(newick_string, newick_format)

    leaves, internal_nodes = node_type(T)

    vector_branches, vector_distances = branch_lengths(T)

    # """get sequences from msa (in phylip format or fasta format) save them into new INPUTMSA.txt file """

    modify_seqfile(pathDATA, leaves, option)

    # """" no idea, where I need this, copy .state to memory.csv"""
    nodes_number, nucleotides_sites = number_nodes_sites(pathDATA)

    if number_rates > 1:
        numbersitesperrate = subsequences(T, pathDATA, epsilon, number_rates, option)
        rates = save_rates(pathDATA, number_rates)
    else:
        rates = 1

    clades1, clades2 = clades(T, t, newick_format, internal_nodes, leaves)

    save_clades(pathDATA, number_rates, clades1, clades2, newick_format, rates)

    if number_rates > 1:
        sequences_clades(
            pathDATA,
            number_rates,
            nodes_number,
            nucleotides_sites,
            clades1,
            clades2,
            option,
            newick_format,
            internal_nodes,
            numbersitesperrate,
        )
    else:
        sequences_clades(
            pathDATA,
            number_rates,
            nodes_number,
            nucleotides_sites,
            clades1,
            clades2,
            option,
            newick_format,
            internal_nodes,
            [],
        )

    (
        state_frequencies_vect,
        model_and_frequency,
    ) = parse_rate_and_frequencies_and_create_model_files(
        pathDATA, number_rates, dimension, model
    )

    """ get the right eigenvector(s) of the dominate non-zero eigenvalue"""
    array_eigenvectors, multiplicity = spectral_decomposition(dimension, pathDATA)

    """ calculate the posterior probabilities using IQ-TREE, see .state files"""
    #run_iqtree_for_each_clade_parallel(pathFOLDER, number_rates, chosen_rate, pathIQTREE)
    run_iqtree_for_each_clade(pathFOLDER, number_rates, chosen_rate, pathIQTREE)

    print(
        "{:6s}\t{:6s}\t{:6s}\t{:6s}\t{:14s}\t{:14s}\t{:100s}".format(
            "Order", "delta", "c_s", "p-value", "Branch status", "T2T status", "Branch"
        )
    )

    """ Calculations of the saturation test for each branch"""
    for i in range(0, len(internal_nodes) + len(leaves)):
       
        """ preparation for results file """
        if i == 0:
            T = Tree(t, format=newick_format)

            results_file = open(
                pathFOLDER + "/resultsRate" + chosen_rate + ".txt", "w"
            ) # To store test results. We open file in first iteration (branch).

            results_file.write("\n")
            results_file.write(
                "{:6s}\t{:6s}\t{:6s}\t{:6s}\t{:14s}\t{:14s}\t{:100s}".format(
                    "Order",
                    "delta",
                    "c_s",
                    "p-value",
                    "Branch status ",
                    "T2T status",
                    "Branch",
                )
            )
            results_file.write("\n")

        """ get the posterior probabilities """
        if number_rates == 1:  # if not gamma model
            # read the posterior probabilities of the left subtree from .state file
            posterior_probabilities_left_subtree = parse_output_state_frequencies(
                f"{pathFOLDER}clades/Branch{str(i)}_clade1/output.state"
            )
            # read the posterior probabilities of the right subtree from .state file
            posterior_probabilities_right_subtree = parse_output_state_frequencies(
                f"{pathFOLDER}clades/Branch{str(i)}_clade2/output.state"
            )

        else:  # if gamma model
            # read the posterior probabilities of the left subtree from .state file
            posterior_probabilities_left_subtree = parse_output_state_frequencies(
                f"{pathFOLDER}/subsequences/subseq{chosen_rate}/clades/Branch{i}_clade1/output.state"
            )
            # read the posterior probabilities of the right subtree from .state file
            posterior_probabilities_right_subtree = parse_output_state_frequencies(
                f"{pathFOLDER}/subsequences/subseq{chosen_rate}/clades/Branch{i}_clade2/output.state"
            )
        
        """ calculation of the test-statistic"""
        if i < len(internal_nodes):
            branch_type = "internal"

        delta, c_s,c_sTwoSequence, p_value, result_test, result_test_tip2tip = calculate_test_statistic(
            multiplicity,
            array_eigenvectors,
            posterior_probabilities_left_subtree,
            posterior_probabilities_right_subtree,
            dimension,
            branch_type,
            alpha,
        )

        results_file.write(
            "{:6d}\t{:6.4f}\t{:6.4f}\t{:6.10f}\t{:14s}\t{:14s}\t{:100s}".format(
                i + 1,
                delta,
                c_s,
                p_value,
                result_test,
                result_test_tip2tip,
                vector_branches[i],
            )
        )
        results_file.write("\n\n")

        print(
            "{:6d}\t{:6.4f}\t{:6.4f}\t{:6.10f}\t{:14s}\t{:14s}\t{:100s}".format(
                i + 1,
                delta,
                c_s,
                p_value,
                result_test,
                result_test_tip2tip,
                vector_branches[i],
            )
        )

    results_file.write(
        "\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is  {:6.4f}".format(
            c_sTwoSequence
        )
    )
    results_file.write(
        "\n\nFor better reference, this is the reconstructed tree topology :\n\n"
    )
    results_file.write(
        T.copy("newick").get_ascii(attributes=["name", "label", "distance"])
    )
    #results_file.close()

    print(
        "\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is ",
        "{:.4f}".format(c_sTwoSequence),
    )

    # TODO New Function but it has to be tested and checked
    # Write results to file and append saturation information to newick string
    write_results_and_newick_tree(
        results_list=results_list,
        newick_string=newick_string,
        path_folder=pathFOLDER,
        chosen_rate=chosen_rate,
        c_sTwoSequence=c_sTwoSequence,
        T=T,
    )


# ========= Mapping Code to Newick String ===========
def map_values_to_newick(value_rows, newick_string):
    # Parsing the file of values into a dictionary
    values_dict = {}

    """
     {
                "index": i + 1,
                "delta": delta,
                "c_s": c_s,
                "p_value": p_value,
                "result_test": result_test,
                "result_test_tip2tip": result_test_tip2tip,
                "vector_branches": vector_branches[i],
    }
    """
    for row in value_rows:
        node_one, node_two = row["vector_branches"].split("-")
        node_two = node_two.split("*")[0]  # remove trailing '*'
        values_dict[node_two] = {
            "delta": row["delta"],
            "c_s": row["c_s"],
            "p-value": row["p_value"],
            "result_test": row["result_test"],
            "status": row["result_test_tip2tip"],
        }

    saturised_newick_string = map_values_to_newick_regex(values_dict, newick_string)

    return saturised_newick_string


def map_values_to_newick_regex(values_dict, newick_string):
    for node_name, values in values_dict.items():
        delta = values["delta"]
        c_s = values["c_s"]
        p_value = values["p-value"]
        branch_status = values["status"]
        newick_string = re.sub(
            rf"({node_name})",
            rf"\1[delta={delta}; c_s={c_s}; p_value={p_value}; branch_status={branch_status}]",
            newick_string,
        )

    return newick_string
