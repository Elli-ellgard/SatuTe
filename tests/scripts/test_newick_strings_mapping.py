from satute.ostream import insert_metadata_into_newick


def test_insert_metadata_with_brackets():
    newick = "(A:0.1[pA=0.2,pB=0.3],B:0.2[pC=0.4,pD=0.5]):0.3;"
    target_node = "A"
    meta_data = "pE=0.6,pF=0.7"
    expected_newick = "(A:0.1[pA=0.2,pB=0.3,pE=0.6,pF=0.7],B:0.2[pC=0.4,pD=0.5]):0.3;"
    assert insert_metadata_into_newick(newick, target_node, meta_data) == expected_newick
    
def test_insert_metadata_multiple_nodes():
    newick = "(A:0.1[pA=0.2,pB=0.3],B:0.2[pC=0.4,pD=0.5],C:0.3[pE=0.6,pF=0.7]):0.3;"
    target_node = "B"
    meta_data = "pG=0.8,pH=0.9"
    expected_newick = "(A:0.1[pA=0.2,pB=0.3],B:0.2[pC=0.4,pD=0.5,pG=0.8,pH=0.9],C:0.3[pE=0.6,pF=0.7]):0.3;"
    assert insert_metadata_into_newick(newick, target_node, meta_data) == expected_newick

def test_insert_metadata_empty_newick():
    newick = ""
    target_node = "A"
    meta_data = "pE=0.6,pF=0.7"
    expected_newick = ""
    assert insert_metadata_into_newick(newick, target_node, meta_data) == expected_newick

def test_insert_metadata_target_node_not_found():
    newick = "(A:0.1[pA=0.2,pB=0.3],B:0.2[pC=0.4,pD=0.5]):0.3;"
    target_node = "D"
    meta_data = "pE=0.6,pF=0.7"
    expected_newick = "(A:0.1[pA=0.2,pB=0.3],B:0.2[pC=0.4,pD=0.5]):0.3;"
    assert insert_metadata_into_newick(newick, target_node, meta_data) == expected_newick