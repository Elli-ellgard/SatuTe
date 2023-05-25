

import numpy as np
import regex as re
import pandas as pd
import ete3
from ete3 import Tree
import csv
import os
import scipy
import scipy.linalg
import subprocess

def remove_filename(path):
    parts = path.split('/')
    nameFILE = parts[-1]
    parts.pop()
    pathFOLDER = '/'.join(parts)+"/" #Path where the we create folders
    return pathFOLDER

def read_tree(path,newickformat):

    filne = path +".iqtree"
    with open(filne, 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "Tree in newick format:" in line:
                t = lines[i+2]
                break

    T = Tree(t,format=newickformat)

    for node in T.traverse("levelorder"):
        l = len(node.name)
        for i in range(len(t)):
            if t[i:i+l]==str(node.name) and (t[i+l]==";" or t[i+l]==":"):
                t = t[:i+l] + '*' + t[i+l:]

    T = Tree(t,format=newickformat)

    root = T.get_tree_root()
    root_children = root.get_children()
    leaves = T.get_leaves()

    root_children_leaves = []
    for i in root_children:
        if i in leaves:
            root_children_leaves.append(i)

    if len(root_children_leaves)>=2:
        for i in root_children_leaves:
            for j in range(len(t)):
                if t[j:j+len(i.name)]==i.name:
                    cont = 0
                    split0 = t.split(i.name,1)[0]
                    split1 = t.split(i.name,1)[1]
                    for k in range(len(split0)):
                        if split0[k]=="(":
                            cont+=1
                    if cont>1:
                        t = t.replace(i.name,"")
                        if (len(split1[0:split1.find(")")])<len(split1[0:split1.find(",")])):
                            t = t.replace(split1[0:split1.find(")")],"")
                            t = t[0] + str(i.name) + split1[0:split1.find(")")] + "," + t[1:-1]
                        else:
                            t = t.replace(split1[0:split1.find(",")],"")
                            t = t[0] + str(i.name) + split1[0:split1.find(",")] + "," + t[1:-1]

    for i in range(len(t)):
        if t[i:i+2]==",)":
            t = t.replace(t[i:i+2],")")

    #print(t)
    T = Tree(t,format=newickformat)
    #print("This is the reconstructed tree topology :\n",T.copy("newick").get_ascii(attributes=["name","label","distance"]))
    return t,T

"""## INTERNAL NODES AND LEAVES"""

def node_type(T):

    leaves = []
    for i in T.get_leaves():
        leaves.append(i.name)

    internal_nodes = []
    for node in T.traverse("levelorder"):
        if (node not in T.get_leaves()) and (node in T.get_tree_root()):
            internal_nodes.append(node.name)

    return leaves,internal_nodes

"""## MODIFYING SEQUENCE FILE

#### We modified the nodes names in order to have indistinguishable names. We need to modify as well the sequence file.
"""

def modify_seqfile(path,leaves,option):

    filesequence = path

    if option==1:

        with open(filesequence, 'r+') as f:
            with open(path+".txt", 'w') as writer:
                lines = f.readlines()
                writer.write(lines[0])
                for l in lines:
                    for i in leaves:
                        if l[0:len(i)-1]==i[0:len(i)-1]:
                            lmod = l[:len(i)-1] + '*' + l[len(i)-1:]
                            writer.write(lmod)

    elif option==2:

        with open(filesequence, 'r+') as f:
            with open(path+".txt", 'w') as writer:
                lines = f.readlines()
                for line in lines:
                    if line[0]==">":
                        for i in leaves:
                            if line[1:-1]==i[0:len(i)-1]:
                                writer.write(">"+i+"\n")
                    else:
                        writer.write(line)

    else:
        print("FILE FORMAT OPTION NOT VALID")

"""## COMPUTING BRANCHES LENGTHS """

def branches_lengths(T):

    vector_branches = []
    vector_distances = []

    internal_nodes = []
    for node in T.traverse("levelorder"):
        if (node not in T.get_leaves()) and (node in T.get_tree_root()):
            internal_nodes.append(node.name)

    for node in T.traverse("levelorder"): #First internal branches.
        children = node.get_children()

        for child in children:
            if child.name in internal_nodes:
                vector_branches.append(node.name+"-"+child.name)
                vector_distances.append(T.get_distance(node,child))

    for node in T.traverse("levelorder"): #Same for external branches.
        children = node.get_children()

        for child in children:
            if child.name not in internal_nodes:
                vector_branches.append(node.name+"-"+child.name)
                vector_distances.append(T.get_distance(node,child))

    return vector_branches,vector_distances

"""## NUMBER OF NODES AND SITES"""

def number_nodes_sites(path):
    filne = path+".state"
    pathFolder = remove_filename(path)
    with open(filne, 'r+') as f:
        with open(pathFolder+'memory.csv', 'w') as writer:
            lines = f.readlines()
            out = lines[8:-1]
            for i in range(len(lines[8:-1])):
                writer.write(lines[i+8])

    df = pd.read_csv(pathFolder+"memory.csv",sep='\t',engine='python')
    file = open(path+".iqtree", "rt")

    nodes_number = len(df['Node'].unique())
    nucleotides_sites = len(df['Site'].unique())

    return nodes_number,nucleotides_sites

"""## MEMORY VECTOR"""

def memory(path):

    filne = path+".state"
    pathFolder = remove_filename(path)
    with open(filne, 'r+') as f:
        with open(pathFolder+'memory.csv', 'w') as writer:
            lines = f.readlines()
            out = lines[8:-1]
            for i in range(len(lines[8:-1])):
                writer.write(lines[i+8])

    df = pd.read_csv(pathFolder+"memory.csv",sep='\t',engine='python')
    file = open(pathFolder+".iqtree", "rt")

    nodes_order = df['Node'].unique()
    nodes_number = len(df['Node'].unique())
    nucleotides_sites = len(df['Site'].unique())

    vector_pi = [] #vector of stationary state nucleotide frequencies
    data = file.readlines()
    for line in data:
        if 'pi(A)' in line:
            vector_pi.append(float(line[10:-1]))
        if 'pi(C)' in line:
            vector_pi.append(float(line[10:-1]))
        if 'pi(G)' in line:
            vector_pi.append(float(line[10:-1]))
        if 'pi(T)' in line:
            vector_pi.append(float(line[10:-1]))

    matrix_norm = [] #matrix to store the norm of each node
    matrix_memory = [] #matrix to store the memory of each node

    for i in range(int(nodes_number)):
        norm_vect = [] #vector (of length=nucleotides_sites) to store for each site in node: {sum(1 to 4 nucleotides)[state(i)^2/stationary_dist(i)]-1}

        for j in range(int(nucleotides_sites)):
            index = i*int(nucleotides_sites)+j-1

            suma = 0
            for k in range(4): #sum in all nucleotides of the current site
                suma += df.iloc[index,k+3]*df.iloc[index,k+3]/vector_pi[k]

            norm_vect.append(suma-1) #adding sum-1 to vector

        suma1 = 0
        suma2 = 0
        for j in norm_vect: #summing all elements in vector -> we get 1 number per node, which is the sum of all numbers of the sites
            suma2 += np.sqrt(j)

        matrix_memory.append(suma1/nucleotides_sites) #dividing sum by total number of nucleotides
        matrix_norm.append(suma2/nucleotides_sites)

    matrix_memory = np.array(matrix_memory)
    matrix_norm = np.array(matrix_norm)

    minimum_memory = nodes_order[np.where(matrix_memory == matrix_memory.min())]
    print("Node with minimum memory: ",minimum_memory[0])

    return minimum_memory

"""## SEPARATING CLADES"""

def clades(T,t,newickformat,internal_nodes,leaves):

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
            if t1[j:j+length]==i:
                for k in reversed(range(len(t1[0:j]))):
                    if t1[k]==")":
                        cont += 1
                    elif t1[k]=="(":
                        cont -= 1
                        aux = k

                    if cont==0:
                        clade1 = t1[k:j+length]+";"
                        clades1.append(clade1)
                        break

        root_clade2 = parent_nodes[internal_nodes.index(i)] #parent node of current node that will be the root of clade 2

        if root_clade2 == root.name:

            clade1 = clade1.replace(";","")
            clade2 = t1.replace(clade1,"")
            for j in range(len(clade2)):
                if clade2[j:j+2]==",:":
                    split = clade2.split(",:",1)[1]
                    if split.find(",(")>0:
                        clade2 = clade2.replace(clade2[j:j+split.find(",(")],'',1)
                    else:
                        clade2 = clade2.replace(clade2[j:j+split.find(")")],'',1)

        else:
            r = root_clade2

            #FIRST WE NEED TO SAVE THE NODES FROM ROOT TO END AND THEIR CORRESPONDING CHILDREN

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
                            split = t1.split(children[k].name + ":",1)[1]
                            if split[0:split.find(")")] < split[0:split.find(",")]:
                                auxiliar.append(children[k].name + ":" + split[0:split.find(")")])
                            else:
                                auxiliar.append(children[k].name + ":" + split[0:split.find(",")])
                            auxiliar_parent.append(r)
                        else:
                            node = T.search_nodes(name=children[k].name)[0]
                            node = node.write(format=newickformat)
                            node = node.replace(";","")
                            auxiliar.append(node)
                            auxiliar_parent.append(r)

                r1 = r
                if r1 != root.name:
                    r = T.search_nodes(name=r)[0].up.name
                    current = r1

            #NOW WE SAVE DISTANCES

            distances_fromroottoend = [""]
            for j in range(1,len(vect_fromroottoend)):
                split = t1.split(vect_fromroottoend[j-1]+":",1)[1]
                if split[0:split.find(")")] < split[0:split.find(",")]:
                    distances_fromroottoend.append(split[0:split.find(")")])
                else:
                    distances_fromroottoend.append(split[0:split.find(",")])

            #LET'S CONSTRUCT THE NEW CLADE

            for j in reversed(vect_fromroottoend):
                index = vect_fromroottoend.index(j)
                if j != root_clade2:
                    if auxiliar_parent.count(j)==2:
                        attach = []
                        for k in range(len(auxiliar_parent)):
                            if auxiliar_parent[k] == j:
                                 attach.append(auxiliar[k])
                        clade2 = "(" + attach[0] + "," + attach[1] + ")" + j + ":" + distances_fromroottoend[index]
                    else:
                        clade2 = "(" + auxiliar[index] + "," + clade2 + ")" + j + ":" + distances_fromroottoend[index]
                else:
                    clade2 = "(" + auxiliar[index] + "," + clade2 + ")" + root_clade2 + ";"

        clades2.append(clade2)

    parent_leaves = []
    for i in leaves:
        node = T.search_nodes(name=i)[0]
        parent = node.up
        parent_leaves.append(parent.name)

    for i in range(len(leaves)):
        clade1 = "(copy_" + leaves[i] + ":0.0001000000," + leaves[i] + ":0.0001000000)ROOT;"
        clades1.append(clade1)

        root_clade2 = parent_leaves[i] #parent node of leaf that will be the root of clade 2

        if root_clade2 == root.name:
            split = t.split(leaves[i],1)[1]
            clade2 = t.replace(leaves[i],"")
            clade2 = clade2.replace(split[0:split.find(",")+1],"",1)

        else:
            r = root_clade2

            #FIRST WE NEED TO SAVE THE NODES FROM ROOT TO END AND THEIR CORRESPONDING CHILDREN

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
                            split = t1.split(children[k].name + ":",1)[1]
                            if split[0:split.find(")")] < split[0:split.find(",")]:
                                auxiliar.append(children[k].name + ":" + split[0:split.find(")")])
                            else:
                                auxiliar.append(children[k].name + ":" + split[0:split.find(",")])
                            auxiliar_parent.append(r)
                        else:
                            node = T.search_nodes(name=children[k].name)[0]
                            node = node.write(format=newickformat)
                            node = node.replace(";","")
                            auxiliar.append(node)
                            auxiliar_parent.append(r)

                r1 = r
                if r1 != root.name:
                    r = T.search_nodes(name=r)[0].up.name
                    current = r1

            #NOW WE SAVE DISTANCES

            distances_fromroottoend = [""]
            for j in range(1,len(vect_fromroottoend)):
                split = t1.split(vect_fromroottoend[j-1]+":",1)[1]
                if split[0:split.find(")")] < split[0:split.find(",")]:
                    distances_fromroottoend.append(split[0:split.find(")")])
                else:
                    distances_fromroottoend.append(split[0:split.find(",")])

            #LET'S CONSTRUCT THE NEW CLADE

            for j in reversed(vect_fromroottoend):
                index = vect_fromroottoend.index(j)
                if j != root_clade2:
                    if auxiliar_parent.count(j)==2:
                        attach = []
                        for k in range(len(auxiliar_parent)):
                            if auxiliar_parent[k] == j:
                                 attach.append(auxiliar[k])
                        clade2 = "(" + attach[0] + "," + attach[1] + ")" + j + ":" + distances_fromroottoend[index]
                    else:
                        clade2 = "(" + auxiliar[index] + "," + clade2 + ")" + j + ":" + distances_fromroottoend[index]
                else:
                    clade2 = "(" + auxiliar[index] + "," + clade2 + ")" + root_clade2 + ";"
        clades2.append(clade2)

    return clades1,clades2

"""## DIVIDING SEQUENCE FILE INTO SUBFILES DEPENDING ON WHICH RATE IS MOST PROBABLE"""

def subsequences(T, path,epsilon,number_rates,option):

    filesequence = path+".txt"

    fileprob = path+".siteprob"

    pathFolder = remove_filename(path)
    with open(fileprob, 'r+') as f:
        with open(pathFolder+'prob.csv', 'w') as writer:
            lines = f.readlines()
            for i in range(len(lines)):
                writer.write(lines[i])

    df = pd.read_csv(pathFolder+"prob.csv",sep='\t',engine='python')
    site_rate =[]
    length = len(df.index)
    for i in range(length):
        probs = []
        for j in range(number_rates):
            probs.append(df.iloc[i][j+1])
        max_value = max(probs)
        if max_value < 1/number_rates+epsilon:
            site_rate.append(0)
        else:
            site_rate.append(probs.index(max_value)+1)
    numbersitesperrate = []
    for j in range(number_rates):
        numbersitesperrate.append(site_rate.count(j+1))

    if option==1:
        for i in range(number_rates):
            os.makedirs(pathFolder + "subsequences/subseq" + str(i+1), exist_ok=True)
            fseq = open(pathFolder + "subsequences/subseq" + str(i+1) + "/sequence.txt", "w+")
            with open(filesequence, 'r+') as f:
                lines = f.readlines()
                l = lines[0]
                fseq.write(l[0:l.find(" ")] + " " + str(site_rate.count(i+1)) + "\n")
                for j in range(1, len(lines)):
                    line = lines[j]
                    fseq.write(line.split(" ",1)[0])
                    seq = line.split(" ",1)[1]
                    for k in range(len(seq)):
                        if seq[k] != " ":
                            index = k
                            break
                    seq = seq[index:-1]
                    #cassius
                    fseq.write((k+1)*" ")
                    for k in range(len(site_rate)):
                        if (site_rate[k]==i+1):
                            fseq.write(seq[k])
                    fseq.write("\n")
    elif option==2:
        leaves = []
        for i in T.get_leaves():
            leaves.append(i.name)
        for i in range(number_rates+1):
            os.makedirs(pathFolder + "subsequences/subseq" + str(i+1), exist_ok=True)
            fseq = open(pathFolder + "subsequences/subseq" + str(i+1) + "/sequence.txt", "w+")
            with open(filesequence, 'r+') as f:
                lines = f.readlines()
                seq = ""
                for j in range(len(lines)):
                    if lines[j][1:-1] in leaves:
                        if(len(seq)>0):
                            for k in range(len(site_rate)):
                                if (site_rate[k]==i+1):
                                    fseq.write(seq[k])
                            fseq.write("\n")
                        fseq.write(lines[j])
                        seq = ""
                    else:
                        seq = seq+lines[j][0:lines[j].find(" ")]
            for k in range(len(site_rate)):
                if (site_rate[k]==i+1):
                    fseq.write(seq[k])
    else:
        print("FILE FORMAT OPTION NOT VALID")

    return numbersitesperrate

"""## SAVING MUTATION RATES (in case of the Gamma model)"""

def save_rates(path,number_rates):

    rates = []
    filne = path+".iqtree"
    with open(filne, 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if " Category  Relative_rate  Proportion" in line:
                if lines[i+1][2:3]==0:
                    for j in range(number_rates):
                        rates.append(float(lines[i+j+2][12:20]))
                else:
                    for j in range(number_rates):
                        rates.append(float(lines[i+j+1][12:20]))

    return rates

"""## MODIFY BRANCHES LENGTHS (if it corresponds), ADD FOO AND SAVE CLADES"""

def save_clades(path,number_rates,clades1,clades2,newickformat,rates):
    pathFolder = remove_filename(path)
    if number_rates==1:
        for i in range(len(clades1)):
            os.makedirs(pathFolder + "clades/Branch" + str(i) + "_clade1/",exist_ok=True)
            os.makedirs(pathFolder + "clades/Branch" + str(i) + "_clade2/",exist_ok=True)

            clades1[i] = "(FOO:0.00000000010," + clades1[i][0:clades1[i].rfind(")")] + ")" + clades1[i][clades1[i].rfind(")"):-1] + ";"
            clades2[i] = "(FOO:0.00000000010," + clades2[i][0:clades2[i].rfind(")")] + ")" + clades2[i][clades2[i].rfind(")"):-1] + ";"

            f1 = open(pathFolder + "clades/Branch" + str(i) + "_clade1/tree.txt", "w")
            f1.write(clades1[i])

            f2 = open(pathFolder + "clades/Branch" + str(i) + "_clade2/tree.txt", "w")
            f2.write(clades2[i])

    else:
        for j in range(number_rates):
            for i in range(len(clades1)):
                cl1 = clades1[i]
                cl2 = clades2[i]

                os.makedirs(pathFolder + "subsequences/subseq" + str(j+1) + "/clades/Branch" + str(i) + "_clade1/",exist_ok=True)
                os.makedirs(pathFolder + "subsequences/subseq" + str(j+1) + "/clades/Branch" + str(i) + "_clade2/",exist_ok=True)

                C1 = Tree(cl1,format=newickformat)
                r = C1.get_tree_root()
                for node in C1.traverse("levelorder"):
                    node.dist = rates[j]*node.dist
                c = C1.write(format=newickformat)
                cl1 = c[0:-1] + r.name + ";"

                C2 = Tree(cl2,format=newickformat)
                r = C2.get_tree_root()
                for node in C2.traverse("levelorder"):
                    node.dist = rates[j]*node.dist
                c = C2.write(format=newickformat)
                cl2 = c[0:-1] + r.name + ";"

                cl1 = "(FOO:0.00000000010," + cl1[0:cl1.rfind(")")] + ")" + cl1[cl1.rfind(")"):-1] + ";"
                cl2 = "(FOO:0.00000000010," + cl2[0:cl2.rfind(")")] + ")" + cl2[cl2.rfind(")"):-1] + ";"

                f1 = open(pathFolder + "subsequences/subseq" + str(j+1) + "/clades/Branch" + str(i) + "_clade1/tree.txt", "w")
                f1.write(cl1)

                f2 = open(pathFolder + "subsequences/subseq" + str(j+1) + "/clades/Branch" + str(i) + "_clade2/tree.txt", "w")
                f2.write(cl2)

"""## MODIFYING BRANCHES LENGTHS IN FULL TREES"""

def modify_fulltree(path,T,rates,newickformat):

    for j in range(rates):
        r = T.get_tree_root()
        for node in T.traverse("levelorder"):
            node.dist = rates[j]*node.dist
        t = T.write(format=newickformat)
        t = t[0:-1] + r.name + ";"

        f = open(path + "subsequences/subseq" + str(j+1) + "/tree.txt", "w")
        f.write(t)

"""## CREATING FOR EACH CLADE THE NUCLEOTIDE SEQUENCE FILE"""

def sequences_clades(path,number_rates,nodes_number,nucleotides_sites,clades1,clades2,option,newickformat,internal_nodes,numbersitesperrate):
    pathFolder = remove_filename(path)

    if number_rates==1:
        if option==1:
            filesequence = path+".txt"

            for i in range(len(clades1)):
                C1 = Tree(clades1[i],format=newickformat)

                f1 = open(pathFolder + "clades/Branch" + str(i) + "_clade1/sequence.txt", "w")

                leaves = []
                for k in C1.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, 'r+') as f:
                    lines = f.readlines()
                    f1.write(str(len(leaves)) + lines[0][lines[0].find(" "):-1])
                    f1.write("\n")
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[0:line.find(" ")] in leaves:
                            f1.write(line)

                if i>=len(internal_nodes): #if clade is only a leaf, we need to rewrite it in the sequence file
                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0:line.find(" ")] in leaves:
                                f1.write("copy_"+line)

                f1.write("FOO  " + "N"*(int(lines[0][lines[0].find(" ")+1:-1])-1) + "A")

            for i in range(len(clades2)):
                C2 = Tree(clades2[i],format=newickformat)

                f2 = open(pathFolder + "clades/Branch" + str(i) + "_clade2/sequence.txt", "w")

                leaves = []
                for k in C2.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, 'r+') as f:
                    lines = f.readlines()
                    f2.write(str(len(leaves)) + lines[0][lines[0].find(" "):-1])
                    f2.write("\n")
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[0:line.find(" ")] in leaves:
                            f2.write(line)

                f2.write("FOO  " + "N"*(int(lines[0][lines[0].find(" ")+1:-1])-1) + "A")

        elif option==2:
            filesequence = path+".txt"

            for i in range(len(clades1)):
                C1 = Tree(clades1[i],format=newickformat)

                f1 = open(pathFolder + "clades/Branch" + str(i) + "_clade1/sequence.txt", "w")

                leaves = []
                for k in C1.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, 'r+') as f:
                    lines = f.readlines()
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[1:-1] in leaves:
                            f1.write(line)
                            for k in range(j+1,len(lines)):
                                line = lines[k]
                                if line[0]==">":
                                    break
                                else:
                                    f1.write(line)

                if i>=len(internal_nodes): #if clade is only a leaf, we need to rewrite it in the sequence file
                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f1.write(line[0]+"copy_"+line[1:-1]+"\n")
                                for k in range(j+1,len(lines)):
                                    line = lines[k]
                                    if line[0]==">":
                                        break
                                    else:
                                        f1.write(line)

                f1.write(">FOO\n" + "N"*(nucleotides_sites-1) + "A")

            for i in range(len(clades2)):
                C2 = Tree(clades2[i],format=newickformat)

                f2 = open(pathFolder + "clades/Branch" + str(i) + "_clade2/sequence.txt", "w")

                leaves = []
                for k in C2.get_leaves():
                    leaves.append(k.name)

                with open(filesequence, 'r+') as f:
                    lines = f.readlines()
                    for j in range(0, len(lines)):
                        line = lines[j]
                        if line[1:-1] in leaves:
                            f2.write(line)
                            for k in range(j+1,len(lines)):
                                line = lines[k]
                                if line[0]==">":
                                    break
                                else:
                                    f2.write(line)

                f2.write(">FOO\n" + "N"*(nucleotides_sites-1) + "A")

        else:
            print("FILE FORMAT OPTION NOT VALID")

    else:
        if option==1:
            for r in range(number_rates):

                filesequence = pathFolder + "subsequences/subseq" + str(r+1) + "/sequence.txt"

                for i in range(len(clades1)):
                    C1 = Tree(clades1[i],format=newickformat)

                    f1 = open(pathFolder + "subsequences/subseq" + str(r+1) + "/clades/Branch" + str(i) + "_clade1/sequence.txt", "w")

                    leaves = []
                    for k in C1.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        f1.write(str(len(leaves)+1) + lines[0][lines[0].find(" "):-1])
                        f1.write("\n")
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0:line.find(" ")] in leaves:
                                f1.write(line)

                    if i>=len(internal_nodes): #if clade is only a leaf, we need to rewrite it in the sequence file
                        with open(filesequence, 'r+') as f:
                            lines = f.readlines()
                            for j in range(0, len(lines)):
                                line = lines[j]
                                if line[0:line.find(" ")] in leaves:
                                    f1.write("copy_"+line)

                    f1.write("FOO  " + "N"*(numbersitesperrate[r]-1) + "A")

                for i in range(len(clades2)):
                    C2 = Tree(clades2[i],format=newickformat)

                    f2 = open(pathFolder + "subsequences/subseq" + str(r+1) + "/clades/Branch" + str(i) + "_clade2/sequence.txt", "w")

                    leaves = []
                    for k in C2.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        f2.write(str(len(leaves)+1) + lines[0][lines[0].find(" "):-1])
                        f2.write("\n")
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[0:line.find(" ")] in leaves:
                                f2.write(line)

                    f2.write("FOO  " + "N"*(numbersitesperrate[r]-1) + "A")

        elif option==2:
            for r in range(number_rates):

                filesequence =  pathFolder + "subsequences/subseq" + str(r+1) + "/sequence.txt"

                for i in range(len(clades1)):
                    C1 = Tree(clades1[i],format=newickformat)

                    f1 = open(pathFolder + "subsequences/subseq" + str(r+1) + "/clades/Branch" + str(i) + "_clade1/sequence.txt", "w")

                    leaves = []
                    for k in C1.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f1.write(line)
                                for k in range(j+1,len(lines)):
                                    line = lines[k]
                                    if line[0]==">":
                                        break
                                    else:
                                        f1.write(line)
                                f1.write("\n")

                    if i>=len(internal_nodes): #if clade is only a leaf, we need to rewrite it in the sequence file
                        with open(filesequence, 'r+') as f:
                            lines = f.readlines()
                            for j in range(0, len(lines)):
                                line = lines[j]
                                if line[1:-1] in leaves:
                                    f1.write(line[0]+"copy_"+line[1:-1]+"\n")
                                    for k in range(j+1,len(lines)):
                                        line = lines[k]
                                        if line[0]==">":
                                            break
                                        else:
                                            f1.write(line)
                                    f1.write("\n")

                    f1.write(">FOO\n" + "N"*(numbersitesperrate[r]-1) + "A")

                for i in range(len(clades2)):
                    C2 = Tree(clades2[i],format=newickformat)

                    f2 = open(pathFolder + "subsequences/subseq" + str(r+1) + "/clades/Branch" + str(i) + "_clade2/sequence.txt", "w")

                    leaves = []
                    for k in C2.get_leaves():
                        leaves.append(k.name)

                    with open(filesequence, 'r+') as f:
                        lines = f.readlines()
                        for j in range(0, len(lines)):
                            line = lines[j]
                            if line[1:-1] in leaves:
                                f2.write(line)
                                for k in range(j+1,len(lines)):
                                    line = lines[k]
                                    if line[0]==">":
                                        break
                                    else:
                                        f2.write(line)
                                f2.write("\n")

                    f2.write(">FOO\n" + "N"*(numbersitesperrate[r]-1) + "A")

        else:
            print("FILE FORMAT OPTION NOT VALID")

"""## RATE PARAMETER AND STATE FREQUENCIES (from the original reconstruction)"""

def rate_and_frequenciesALTERNATIVE(path,number_rates, dimension):

    with open(path+".iqtree", "r") as f:
        found = 0
        number_lines = 0
        modelString = ""
        for line in f:
            if "Rate parameter R:" in line:
                found = 1
            if found and number_lines < dimension*(dimension-1)//2+1:
                if number_lines > 1:
                    modelString += line[7:]
                number_lines+=1
        separated = modelString.splitlines()
        modelFinal = "GTR{"
        modelFinal += separated[0]
        for words in separated[1:]:
            modelFinal += "," + words
        modelFinal += "}"

    with open(path+".iqtree", "r") as f:
        #FIX: Use split instead of character position, too delicate otherwise
        found = 0
        number_lines = 0
        frequencyString = ""
        for line in f:
            if "State frequencies:" in line:
                found = 1
            if found and number_lines < dimension+2:
                if "equal frequencies" in line:
                    frequencyString = "0.25\n"*dimension
                    break
                if number_lines > 1:
                    frequencyString += line[10:]
                number_lines+=1
        separated = frequencyString.splitlines()
        frequencyFinal = "FU{"
        frequencyFinal += separated[0]
        for words in separated[1:]:
            frequencyFinal += "," + words
        frequencyFinal += "}"
        modelAndFrequency = modelFinal+"+"+frequencyFinal
    state_frequencies_vect = []
    for words in separated:
        state_frequencies_vect.append(float(words))
    pathFolder = remove_filename(path)
    if number_rates==1:
        f1 = open(pathFolder + "model.txt", "w")
        f1.write(modelAndFrequency)
    else:
        for i in range(number_rates):
            f1 = open(pathFolder + "subsequences/subseq" + str(i+1) + "/model.txt", "w")
            f1.write(modelAndFrequency)
    return state_frequencies_vect

"""## DIAGONALISATION OF THE RATE MATRIX

"""

def diagonalisation(n,path):

    ratematrix = np.zeros((n,n))
    phimatrix = np.zeros((n,n))

    filne = path+".iqtree"
    with open(filne, 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "Rate matrix Q:" in line:
                for j in range(n):
                    listEntries = lines[i+j+2].split()
                    for k in range(n):
                        if ('e' in listEntries[k+1]):
                            ratematrix[j,k] = "0"
                        else:
                            ratematrix[j,k]= listEntries[k+1]
            if "State frequencies: (empirical counts from alignment)" in line:
                for j in range(n):
                    entry = lines[i+j+2].split()
                    phimatrix[j,j] = entry[2]
            elif "State frequencies: (equal frequencies)" in line:
                for j in range(n):
                    phimatrix[j,j] = 0.25

    M = scipy.linalg.fractional_matrix_power(phimatrix, +1/2)@ratematrix
    M = M@scipy.linalg.fractional_matrix_power(phimatrix, -1/2)

    lamb,w =  np.linalg.eig(M) #Compute the eigenvalues and right eigenvectors.
    idx = lamb.argsort()[::-1] #Order from large to small.
    lamb = lamb[idx]
    w = w[:,idx]

    lamb_nozero = [] #list of eigenvalues without 0
    for i in lamb:
        if i>0.00999 or i<-0.00999:
            lamb_nozero.append(i)

    index = []
    max_lambda = max(lamb_nozero)
    index.append((lamb.tolist()).index(max_lambda))

    lamb_nozero.remove(max_lambda)
    while len(lamb_nozero)>0:
        max_lambda_it = max(lamb_nozero)
        if abs(max_lambda_it-max_lambda)<0.01:
            index.append((lamb.tolist()).index(max_lambda_it))
        lamb_nozero.remove(max_lambda_it)

    array_eigenvectors = []

    v1 = scipy.linalg.fractional_matrix_power(phimatrix, -1/2)@w[:,index[0]]
    h1 = scipy.linalg.fractional_matrix_power(phimatrix, +1/2)@w[:,index[0]]

    array_eigenvectors.append(v1)

    multiplicity = len(index)
    if multiplicity>1:
        for i in range(1,multiplicity):
            v1 = scipy.linalg.fractional_matrix_power(phimatrix, -1/2)@w[:,index[i]]
            h1 = scipy.linalg.fractional_matrix_power(phimatrix, +1/2)@w[:,index[i]]
            array_eigenvectors.append(v1)

    return array_eigenvectors,multiplicity


def saturationTest(pathDATA, pathIQTREE, runIQTREE = True, runBOOTSRAP = True, dimension = 4, number_rates = 4, chosen_rate = str(4), z_alpha = 2.33, newickformat = 1, epsilon = 0.01, rawMemory = True):
    """IQ-TREE"""
    if runIQTREE:
        modelRatesInv = "GTR+G"+str(number_rates)+"+I"
        modelHomo = "GTR"
        if runBOOTSRAP:
            if number_rates>1:
                subprocess.run([pathIQTREE, "-redo", "-s", pathDATA, "-m", modelRatesInv, "-asr", "-wspr", "-seed", "10", "-quiet", "-b", "100"])
            else:
                subprocess.run([pathIQTREE, "-redo", "-s", pathDATA, "-m", modelHomo, "-asr", "-seed", "10", "-quiet", "-b", "100"])

        else:
            if number_rates>1:
                subprocess.run([pathIQTREE, "-redo", "-s", pathDATA, "-m", modelRatesInv, "-asr", "-wspr", "-seed", "10", "-quiet"])
            else:
                subprocess.run([pathIQTREE, "-redo", "-s", pathDATA, "-m", modelHomo, "-asr", "-seed", "10", "-quiet"])

    """Check type if the alignment is Phylip or fasta format"""

    option = 1 #If sequence file has phylip format
    # Check if otherwise we assume fasta file:
    with open(pathDATA, "r") as toData:
        firstLine = toData.readline().strip('\n')
        words = firstLine.split()
        for word in words:
            if not word.isnumeric():
                option = 2 #Since the first line is not only numbers, we assume it is a fasta file

    """PATHS"""
    pathFOLDER = remove_filename(pathDATA)
    "MAIN"
    t,T =  read_tree(pathDATA,newickformat)
    leaves,internal_nodes =  node_type(T)
    vector_branches,vector_distances =  branches_lengths(T)

    modify_seqfile(pathDATA,leaves,option)

    nodes_number,nucleotides_sites =  number_nodes_sites(pathDATA)

    if number_rates>1:
        numbersitesperrate =  subsequences(T,pathDATA,epsilon,number_rates,option)
        rates =  save_rates(pathDATA,number_rates)
    else:
        rates = 1

    clades1,clades2 =  clades(T,t,newickformat,internal_nodes,leaves)
    save_clades(pathDATA,number_rates,clades1,clades2,newickformat,rates)
    if number_rates> 1:
        sequences_clades(pathDATA,number_rates,nodes_number,nucleotides_sites,clades1,clades2,option,newickformat,internal_nodes, numbersitesperrate)
    else:
        sequences_clades(pathDATA,number_rates,nodes_number,nucleotides_sites,clades1,clades2,option,newickformat,internal_nodes, [])
    state_frequencies_vect =  rate_and_frequenciesALTERNATIVE(pathDATA,number_rates, dimension)

    array_eigenvectors,multiplicity =  diagonalisation(dimension,pathDATA)



    """BASH SCRIPT BEFORE TEST"""
    if (number_rates>1):
        pathNEWFOLDER  = pathFOLDER+"subsequences/subseq"+chosen_rate+"/clades/*"

        with open(pathFOLDER+"subsequences/subseq"+chosen_rate+"/model.txt", "r") as toModel:
            modelAndFrequency = toModel.readline().strip('\n')
        script = """
                for d in """ + pathNEWFOLDER + """; do
                    cd "$d"
                    """ +pathIQTREE +  """ -s sequence.txt -te tree.txt -m \"'\"""" + modelAndFrequency+ """\"'\" -asr -blfix -o FOO -pre output -redo -quiet
                done   
                """
    else:
        pathNEWFOLDER  = pathFOLDER+"clades/*"
        with open(pathFOLDER+"model.txt", "r") as toModel:
            modelAndFrequency = toModel.readline().strip('\n')
        script = """
                for d in """ + pathNEWFOLDER + """; do
                    cd "$d"
                    """ +pathIQTREE +  """ -s sequence.txt -te tree.txt -m \"'\"""" + modelAndFrequency+ """\"'\" -asr -blfix -o FOO -pre output -redo -quiet
                done"""
                
    os.system("bash -c '%s'" % script)



    """SATURATION TEST FOR ALL BRANCHES"""

    U = 1.0/float(min(state_frequencies_vect)) - 1
    K = dimension - 1
    number_standard_deviations = 2 #Confidence intervals of 98% (one sided)


    print("{:6s}  {:6s}  {:6s}  {:14s} {:14s} {:100s}".format("Order"," delta"," c_s","Branch status","T2T status", " Branch","\n"))
    for i in range(0,len(internal_nodes)+len(leaves)):

        if number_rates==1: #if not gamma model

            file1 = pathFOLDER+"clades/Branch"+str(i)+"_clade1/output.state"
            with open(file1, 'r+') as f1:
                with open(pathFOLDER+"clades/Branch"+str(i)+"_clade1/memory.csv", 'w') as writer:
                    lines = f1.readlines()
                    out = lines[8:-1]
                    for j in range(len(lines[8:])):
                        writer.write(lines[j+8])

            file2 = pathFOLDER+"clades/Branch"+str(i)+"_clade2/output.state"
            with open(file2, 'r+') as f2:
                with open(pathFOLDER+"clades/Branch"+str(i)+"_clade2/memory.csv", 'w') as writer:
                    lines = f2.readlines()
                    out = lines[8:-1]
                    for j in range(len(lines[8:])):
                        writer.write(lines[j+8])

            df1 = pd.read_csv(pathFOLDER+"clades/Branch"+str(i)+"_clade1/memory.csv",sep='\t',engine='python')
            df2 = pd.read_csv(pathFOLDER+"clades/Branch"+str(i)+"_clade2/memory.csv",sep='\t',engine='python')
            number_sites = len(df1['Site'].unique())
            number_nodes_1 = len(df1['Node'].unique())
            number_nodes_2 = len(df2['Node'].unique())

            if i==0:
                T = Tree(t,format=newickformat)

                results_file = open(pathFOLDER +"/resultsRate"+chosen_rate+".txt", "w") #To store test results. We open file in first iteration (branch).
                #results_file.write(T.copy("newick").get_ascii(attributes=["name","label","distance"]))
                results_file.write("\n")
                results_file.write('{:6s}\t{:6s}\t{:6s}\t{:14s}\t{:14s}\t{:100s}'.format("Order"," delta"," c_s"," Branch status ", "T2T status", "Branch","\n"))
                results_file.write("\n")

        else: #if gamma model

            file1 = pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade1/output.state"
            with open(file1, 'r+') as f1:
                with open(pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade1/memory.csv", 'w') as writer:
                    lines = f1.readlines()
                    out = lines[8:-1]
                    for j in range(len(lines[8:])):
                        writer.write(lines[j+8])

            file2 = pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade2/output.state"
            with open(file2, 'r+') as f2:
                with open(pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade2/memory.csv", 'w') as writer:
                    lines = f2.readlines()
                    out = lines[8:-1]
                    for j in range(len(lines[8:])):
                        writer.write(lines[j+8])

            df1 = pd.read_csv(pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade1/memory.csv",sep='\t',engine='python')
            df2 = pd.read_csv(pathFOLDER + "subsequences/subseq" + chosen_rate + "/clades/Branch" + str(i) + "_clade2/memory.csv",sep='\t',engine='python')
            number_sites = len(df1['Site'].unique())
            number_nodes_1 = len(df1['Node'].unique())
            number_nodes_2 = len(df2['Node'].unique())

            if i==0:
                T = Tree(t,format=newickformat)
                results_file = open(pathFOLDER + "/resultsRate"+chosen_rate+".txt", "w") #To store test results.  We open file in first iteration (branch).
                results_file.write("\n")
                results_file.write('{:6s}\t{:6s}\t{:6s}\t{:14s}\t{:14s}\t{:100s}\n'.format("Order"," delta"," c_s","Branch status","T2T status"," Branch", "\n"))
                results_file.write("\n")

        estimation_dt = np.sqrt(U*min(K,U/4)/number_sites)
        if (not rawMemory):
            upper_ci = number_standard_deviations*estimation_dt
        else:
            upper_ci = float(0)

        if multiplicity == 1: #if D=1

            v1 = array_eigenvectors[0]

            a = [] #vector to store all products v1*rootsitesposteriorprobabilitiescladeA
            b = [] #vector to store all products v1*rootsitesposteriorprobabilitiescladeB

            for k in range(number_sites*(number_nodes_1-1),number_sites*number_nodes_1):
                a.append(v1@np.asarray(df1.iloc[k,3:7]))


            for k in range(number_sites*(number_nodes_2-1),number_sites*number_nodes_2):
                b.append(v1@np.asarray(df2.iloc[k,3:7]))

            delta = np.asarray(a)@np.asarray(b)/number_sites #computing the dominant sample coherence
            #print(np.multiply(a, b))
            if i<len(internal_nodes):
                M_a = np.asarray(a)@np.asarray(a)/number_sites+upper_ci
                M_a = min(1, M_a)
            else: #if clade A is a single leaf
                M_a = 1

            M_b = np.asarray(b)@np.asarray(b)/number_sites+upper_ci
            M_b = min(1, M_b)
            aux = M_a*M_b
            c_s = z_alpha*np.sqrt(aux)/np.sqrt(number_sites) #computing the saturation coherence
            c_sTwoSequence = z_alpha/np.sqrt(number_sites) #computing the saturation coherence between two sequences
        else:
            c_sTwoSequence = multiplicity*z_alpha/np.sqrt(number_sites) #computing the saturation coherence between two sequences
            delta = 0

            for j in range(multiplicity):

                a = []
                b = []

                v1 = array_eigenvectors[j]

                for k in range(number_sites*(number_nodes_1-1),number_sites*number_nodes_1):
                    a.append(v1@np.asarray(df1.iloc[k,3:7]))


                for k in range(number_sites*(number_nodes_2-1),number_sites*number_nodes_2):
                    b.append(v1@np.asarray(df2.iloc[k,3:7]))

                delta += np.asarray(a)@np.asarray(b)

            delta = delta/number_sites

            variance = 0

            for j in range(multiplicity):
                for k in range(multiplicity):

                    a = []
                    b = []

                    v_j = array_eigenvectors[j]
                    v_k = array_eigenvectors[k]

                    for l in range(number_sites*(number_nodes_1-1),number_sites*number_nodes_1):
                        a.append(v_j@np.asarray(df1.iloc[l,3:7]))

                    for l in range(number_sites*(number_nodes_2-1),number_sites*number_nodes_2):
                        b.append(v_k@np.asarray(df2.iloc[l,3:7]))

                    #variance = np.asarray(a)@np.asarray(b)
                    variance += max(np.asarray(a-upper_ci)@np.asarray(b-upper_ci),np.asarray(a+upper_ci)@np.asarray(b+upper_ci),np.asarray(a+upper_ci)@np.asarray(b-upper_ci),np.asarray(a-upper_ci)@np.asarray(b+upper_ci))

            variance = variance/(number_sites*number_sites)

            if variance<0:
                print("VARIANCE ESTIMATION IS NEGATIVE - CONSIDER INCREASING THE NUMBER OF STANDARD DEVIATIONS (number_standard_deviations) (CONFIDENCE INTERVAL)")
                c_s = 999999999
            else:
                c_s = z_alpha*np.sqrt(variance)


        if c_s>delta:
            result_test = "Saturated"
        else:
            result_test = "Informative"

        if c_sTwoSequence>delta:
            result_test_tip2tip = "SatuT2T"
        else:
            result_test_tip2tip = "InfoT2T"
        print("{:6d}  {:6.4f}  {:6.4f}  {:14s} {:14s} {:100s}".format(i+1,delta,c_s,result_test,result_test_tip2tip,vector_branches[i],"\n"))
        results_file.write('{:6d}\t{:6.4f}\t{:6.4f}\t{:14s}\t{:14s}\t{:100s}'.format(i+1,delta,c_s,result_test,result_test_tip2tip,vector_branches[i],"\n"))
        results_file.write("\n\n")

    results_file.write("\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is  {:6.4f}".format(c_sTwoSequence))

    results_file.write("\n\nFor better reference, this is the reconstructed tree topology :\n\n")
    results_file.write(T.copy("newick").get_ascii(attributes=["name","label","distance"]))
    print("\n\nThe T2T status uses as threshold the saturation coherence between two sequences, which is ", "{:.4f}".format(c_sTwoSequence))


if __name__ == '__main__':

    print("Hello! If you give me an alignment, I can check how saturated it is.\n"
          "I will divide the alignment into regions using IQ-TREE to reconstruct the phylogeny assuming the GTR+I+Gamma model.\n"
          "Then I will separate all sites depending on how fast they evolve.\n"
          "For each of these regions, in the reconstructed phylogeny we will test for branch saturation and tip-to-tip (T2T) saturation.\n")

#    pathDATA = input("What is the path to the alignment?\n")
#    pathIQTREE = input("What is the path to IQ-TREE?\n")
#    numberRates = int(input("How many rate categories do you want to use?\n"))
    
    DIR = "/home/elgert/Desktop/Neue_Daten_2023_05_12/results/pbmc3k/satute/"

    pathDATA = DIR + "/subsampled_alignment_pbmc3k_500hvg_genewisenormed/"+"/subsampled_alignment_pbmc3k_500hvg_genewisenormed.fasta"
    #"/home/elgert/Desktop/Neue_Daten_2023_05_12/satute/subsampled_alignment_pbmc3k_500hvg_genewisenormed/fourRates/subsampled_alignment_pbmc3k_500hvg_genewisenormed.fasta"
    # "/home/elgert/Desktop/Neue_Daten_2023_05_12/satute/subsampled_alignment_pbmc3k_100hvg_genewisenormed/subsampled_alignment_pbmc3k_100hvg_genewisenormed.fasta"
    #"/home/elgert/Desktop/Neue_Daten_2023_05_12/satute/subsampled_alignment_pbmc3k_100randomgenes_genewisenormed/subsampled_alignment_pbmc3k_100randomgenes_genewisenormed.fasta"
#    pathTREE= "/home/elgert/Desktop/Neue_Daten_2023_05_12/results/pbmc3k/subsampled_alignment_pbmc3k_100hvg/boot/subsampled_alignment_pbmc3k_100hvg.fasta.treefile"
    pathIQTREE = "/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2"
    numberRates = 1 

    # saturationTest(pathDATA, pathIQTREE, runIQTREE = True, runBOOTSRAP = True, dimension = 4, number_rates = 4, chosen_rate = str(4), z_alpha = 2.33, newickformat = 1, epsilon = 0.01, rawMemory = True)
    for i in range(numberRates):
        print("Here comes the ", i+1, "'th fastest evolving region: ")
        if i==0:
            saturationTest(pathDATA, pathIQTREE, True, False, 4, numberRates, str(numberRates-i), 2.33, 1, 0.01, True)
        else:
            saturationTest(pathDATA, pathIQTREE, False, False, 4, numberRates, str(numberRates-i), 2.33, 1, 0.01, True)

