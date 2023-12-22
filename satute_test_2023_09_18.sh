iqtree=iqtree2
python=python3

##  Test cases from Clemens
#DIR=./test/Clemens/example_1
#msa=ENSG00000087460_GNAS.fasta

#DIR=./test/Clemens/example_2
#msa=sim-JC+G-alpha1.2-taxa64-len1000bp-bla0.01-blb0.2-blc0.1-rep01.phy

#DIR=./test/Clemens/example_3
#msa=sim-JC+G-AC1-AG1-AT1-CG1-CT1-GT1-alpha1.2-taxa64-len1000bp-bla0.01-blb0.8-blc0.2-rep01.fasta

#DIR=./test/Clemens/example_4
#msa=example.txt

#DIR=./test/Clemens/example_sym_1
#msa=ENSG00000119574_ZBTB45.fasta

DIR=./test/Clemens/exampl_sym_3
msa=PF04055.fasta

#DIR=./test/Clemens/toy_example_GTR+G4
#msa=toy_example_ntaxa_7_run_1-alignment.phy

DIR=./test/Clemens/toy_example_JC
msa=toy_example_ntaxa_7_run_5-alignment.phy

PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi

echo " ============= MODI MSA ===================="

echo ""
echo "---------------------------------------------------------"
echo "TEST 1a: only msa, iqtree exists"
echo "---------------------------------------------------------"
$python satute_cli.py -iqtree $iqtree -msa $DIR/$msa  -alpha 0.05
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 1b: only msa, ufboot"
echo "---------------------------------------------------------"
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -iqtree $iqtree -msa $DIR/$msa -ufboot 1000 -alpha 0.05
echo ""


# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 1c: only msa, boot"
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -iqtree $iqtree -msa $DIR/$msa -boot 100 -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 1d: only msa, iqtree exists, path to executable"
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# path_iqtree="/home/elgert/IQ-TREE/build2/iqtree2"
# $python satute_cli.py -iqtree $path_iqtree -msa $DIR/$msa -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 2a: no msa, no dir => error"
# echo "---------------------------------------------------------"
# $python satute_cli.py -iqtree $iqtree -model JC -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 2b: only msa, iqtree not exists => error"
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -iqtree iqtree3 -msa $DIR/$msa -alpha 0.05
# echo ""


echo " ============= MODI MSA+MODEL ===================="

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3a: msa + model, number of rate categories set"
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa  -model JC+G4 -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3b: msa + model with parameter "
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa  -model GTR{4.39/5.30/4.39/1.0/12.1} -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3c: msa + model with parameter + shape parameter of Gamma model"
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa  -model TIM2{4.39/5.30/12.1}+G{0.7} -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3c: msa + model with parameter +  "
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa  -model TIM2{4.39/5.30/12.1}+R2{0.7/0.1/0.3/0.9} -alpha 0.05
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3d: msa + model, option ufboot"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -model JC+G4 -ufboot 1000

# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 3e: msa + model, option boot"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -model JC+G4 -boot 100
# echo ""


# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 4: msa + model + nr "
# echo "---------------------------------------------------------"
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa  -model JC+G -nr 2 -alpha 0.05
# echo ""


echo " ============= MODI MSA+TREE ===================="

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 5a: msa  + tree  => error"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     mv $DIR/${msa}.treefile $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
#     mv $PDIR/${msa}.treefile $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -tree $DIR/${msa}.treefile
# echo ""

echo " ============= MODI MSA+MODEL+TREE ===================="

echo ""
echo "---------------------------------------------------------"
echo "TEST 5: msa + model + tree"
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
#create output
$python satute_cli.py -msa $DIR/$msa
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    mv $DIR/${msa}.treefile $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
    mv $PDIR/${msa}.treefile $DIR
fi
$python satute_cli.py -msa $DIR/$msa -model JC+R5 -tree $DIR/${msa}.treefile
echo ""
exit
# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 5: msa + model + tree, option boot => error"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     mv $DIR/${msa}.treefile $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
#     mv $PDIR/${msa}.treefile $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -model JC+G4 -tree $DIR/${msa}.treefile -boot 100
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 5: msa + model + tree, option ufboot => error"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     mv $DIR/${msa}.treefile $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
#     mv $PDIR/${msa}.treefile $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -model JC+G4 -tree $DIR/${msa}.treefile -ufboot 1000
# echo ""


# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 6: directory does not exist => error"
# echo "---------------------------------------------------------"

# $python satute_cli.py -dir ${DIR}_bla_bla_bla
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 7: directory  exists, iqtree does not exists(not necessary)"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}*satute.tree ]; then
#     rm -r $DIR/${msa}*satute*
# fi
# $python satute_cli.py -dir $DIR -iqtree iqtree3
# echo ""


# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 8: dir  + msa "
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}*satute.tree ]; then
#     rm -r $DIR/${msa}*satute*
# fi
# $python satute_cli.py -dir $DIR -msa $DIR/$msa
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 7: dir missing siteprob "
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $iqtree  -s $DIR/$msa -m JC+R2 -quiet
# $python satute_cli.py -dir $DIR
# echo ""


echo " ============= Option -edge ===================="

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 5: msa + model + tree, option edge"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     mv $DIR/${msa}.treefile $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
#     mv $PDIR/${msa}.treefile $DIR
# fi
# $python satute_cli.py -msa $DIR/$msa -model JC -tree $DIR/${msa}.treefile -edge "(Node2*, Node1*)"
# echo ""

# echo ""
# echo "---------------------------------------------------------"
# echo "TEST 5: dir, option edge"
# echo "---------------------------------------------------------"
# #clean folder
# PDIR=$(dirname $DIR)
# if [ -e $DIR/${msa}.iqtree ]; then
#     mv $DIR/$msa $PDIR
#     rm -r $DIR
#     mkdir $DIR
#     mv $PDIR/$msa $DIR
# fi
# #create output
# $python satute_cli.py -msa $DIR/$msa

# $python satute_cli.py -dir $DIR -edge "(Node2*, Node1*)"
# echo ""

#python3 satute_cli.py -iqtree $iqtree -msa $alignment -tree random_generated_tree.tree -model $model -alpha 0.05

echo " ============= Option -category ===================="

echo ""
echo "---------------------------------------------------------"
echo "TEST 10a: msa + model , option invalid category => error"
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -msa $DIR/$msa -model JC+G4 -category 5
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 10b: msa + model , option valid, but empty category "
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -msa $DIR/$msa -model JC+G4 -category 2
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 10b: msa + model , option valid category "
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -msa $DIR/$msa -model JC+G4 -category 3
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 10c: directory, option invalid category => error"
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
#create output
$python satute_cli.py -msa $DIR/$msa -model GTR
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}*satute.tree ]; then
    rm -r $DIR/${msa}*satute*
fi
$python satute_cli.py -dir $DIR -category 4
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 10c: directory, option valid category"
echo "---------------------------------------------------------"
#clean folder
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
#create output
$python satute_cli.py -msa $DIR/$msa -model GTR+G6
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}*satute.tree ]; then
    rm -r $DIR/${msa}*satute*
fi
$python satute_cli.py -dir $DIR -category 1
echo ""

#cd visualisation
#./script-annotate-trees-from-satute.sh ../$DIR    
#cd ..
