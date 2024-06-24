
from tests.scripts.satute_test_utils import (
    run_satute,    
    create_destination_dir,
    copy_files_to_dest_dir,
    check_satute_files_dir,    
    filter_files_by_last_name,
    list_filenames_in_directory,    
)

from tests.scripts.fixtures import *

def test_3(dir_aa_path, iqtree, python, satute):
    source_path  = dir_aa_path[0] 
    suffix = "Amino Acid: Directory Test Three Without IQ-TREE File"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
    iq_tree_run_file_list = filter_files_by_last_name(iq_tree_run_file_list, 'iqtree')        
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    # Satute run    
    # Capture sys.exit using pytest.raises
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                '-dir',
                str(dest_dir_path.absolute()),
                "-alpha",
                alpha,
                "-asr",
                "-quiet"
            ])
    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure

def test_6(dir_aa_path, iqtree, python, satute):
    source_path  = dir_aa_path[0] 
    suffix = "Amino Acid: Directory Test 6 dir option edge"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    run_satute(
        [
            '-dir',
            str(dest_dir_path.absolute()),
            "-alpha",
            alpha,
            "-asr",
            "-edge", 
            "(Node2*, Node1*)"
        ])
    # check the files   
    assert check_satute_files_dir(dest_dir_path, categories, alpha, asr)

def test_7(dir_category_aa_path, iqtree, python, satute):
    source_path  = dir_category_aa_path[0] 
    suffix = "Amino Acid: Directory Test 7 With Rate Categories And Edge"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = [1,2,3,4]
    alpha = str(0.05)
    asr = True
    
    run_satute(
        [
            '-dir',
            str(dest_dir_path.absolute()),
            "-alpha",
            alpha,
            "-asr",
            "-edge", 
            "(Node2*, Node1*)"
        ])
    # check the files   
    assert check_satute_files_dir(dest_dir_path, categories, alpha, asr)

def test_8(dir_category_aa_path, iqtree, python, satute):
    source_path  = dir_category_aa_path[0] 
    suffix = "Amino Acid: Directory Test 8 With Rate Categories And Edge"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)

        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)
        
    categories = [1,2,3,4]
    alpha = str(0.05)
    asr = True
    
    run_satute(
        [
            '-dir',
            str(dest_dir_path.absolute()),
            "-alpha",
            alpha,
            "-asr",
            "-edge", 
            "(Node2*, Node1*)"
        ])
    # check the files   
    assert check_satute_files_dir(dest_dir_path, categories, alpha, asr)
    
def test_9(dir_category_aa_path, iqtree, python, satute):
    source_path  = dir_category_aa_path[0] 
    suffix = "Amino Acid: Directory Test 9 Without Siteprob file"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
    iq_tree_run_file_list = filter_files_by_last_name(iq_tree_run_file_list, 'siteprob')        
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)
    alpha = str(0.05)

   # Satute run    
    with pytest.raises(SystemExit) as excinfo: 
       run_satute(
           [
               '-dir',
               str(dest_dir_path.absolute()),
               "-alpha",
               alpha,
               "-asr",
           ])

    # Verify the exit code if needed
    assert excinfo.value.code != 0  # assuming exit code 1 for failure    
