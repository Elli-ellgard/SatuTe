DIR=./test/Clemens/toy_example_JC
msa=toy_example_ntaxa_7_run_5-alignment.phy
iqtree=iqtree
python=python3


PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
fi

#create output
$python satute_cli.py -msa $DIR/$msa -iqtree iqtree

PDIR=$(dirname $DIR)
if [ -e $DIR/${msa}.iqtree ]; then
    mv $DIR/$msa $PDIR
    mv $DIR/${msa}.treefile $PDIR
    rm -r $DIR
    mkdir $DIR
    mv $PDIR/$msa $DIR
    mv $PDIR/${msa}.treefile $DIR
fi

$python satute_cli.py -msa $DIR/$msa -model JC+G4 -tree $DIR/${msa}.treefile -iqtree iqtree
echo ""
