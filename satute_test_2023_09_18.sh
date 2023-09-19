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

#DIR=./test/Clemens/example_sym_2
#msa=ENSG00000138316_ADAMTS14.fasta

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

echo ""
echo "---------------------------------------------------------"
echo "TEST 1a: only msa, iqtree exists"
echo "---------------------------------------------------------"
$python satute_cli.py -iqtree $iqtree -msa $DIR/$msa -alpha 0.05
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 1b: only msa, iqtree exists, path to executable"
echo "---------------------------------------------------------"
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
path_iqtree="/home/elgert/IQ-TREE/build2/iqtree2"
$python satute_cli.py -iqtree $path_iqtree -msa $DIR/$msa -alpha 0.05
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 1c: no msa, no dir => error"
echo "---------------------------------------------------------"
$python satute_cli.py -iqtree $iqtree -model JC -alpha 0.05
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 2: only msa, iqtree not exists => error"
echo "---------------------------------------------------------"
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -iqtree iqtree3 -msa $DIR/$msa -alpha 0.05
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 3: msa + model "
echo "---------------------------------------------------------"
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -msa $DIR/$msa  -model GTR+G4 -alpha 0.05
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 4: msa + model + nr "
echo "---------------------------------------------------------"
PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi
$python satute_cli.py -msa $DIR/$msa  -model JC+G -nr 2 -alpha 0.05
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 5a: msa  + tree  => error"
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
$python satute_cli.py -msa $DIR/$msa -tree $DIR/${msa}.treefile
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 5: msa + model"
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
$python satute_cli.py -msa $DIR/$msa -model JC+G4 
echo ""

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
$python satute_cli.py -msa $DIR/$msa -model JC+G4 -tree $DIR/${msa}.treefile
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 6: directory does not exist => error"
echo "---------------------------------------------------------"

$python satute_cli.py -dir ${DIR}_bla_bla_bla
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 7: directory  exists, iqtree does not exists(not necessary)"
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
if [ -e $DIR/${msa}*satute.tree ]; then
    rm -r $DIR/${msa}*satute*
fi
$python satute_cli.py -dir $DIR -iqtree iqtree3
echo ""


echo ""
echo "---------------------------------------------------------"
echo "TEST 8: dir  + msa "
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
if [ -e $DIR/${msa}*satute.tree ]; then
    rm -r $DIR/${msa}*satute*
fi
$python satute_cli.py -dir $DIR -msa $DIR/$msa
echo ""

echo ""
echo "---------------------------------------------------------"
echo "TEST 7: dir missing siteprob "
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
$iqtree  -s $DIR/$msa -m JC+R2 -quiet
$python satute_cli.py -dir $DIR
echo ""






#python3 satute_cli.py -iqtree $iqtree -msa $alignment -tree random_generated_tree.tree -model $model -alpha 0.05


#cd visualisation
#./script-annotate-trees-from-satute.sh ../$DIR    
#cd ..
