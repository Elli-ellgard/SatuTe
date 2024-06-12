import os

from tests.scripts.satute_test_utils import (
    create_destination_dir,
    copy_files_to_dest_dir,
    run_satute,
    list_filenames_in_directory,
    check_satute_files_dir,
    filter_files_by_last_name
)

from tests.scripts.fixtures import *

def test_1(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test One Everything is the folder Without IQ-Tree"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    # Satute run    
    run_satute(
        [
            '-dir',
            str(dest_dir_path.absolute()),
            "-alpha",
            alpha,
            "-asr",
            "-quiet"
        ])

    # check the files   
    assert check_satute_files_dir(dest_dir_path, categories, alpha, asr)
    
def test_2(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test Two Without Tree File"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
    iq_tree_run_file_list = filter_files_by_last_name(iq_tree_run_file_list, 'tree')        
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    print(iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    # Satute run    
    run_satute(
        [
            '-dir',
            str(dest_dir_path.absolute()),
            "-alpha",
            alpha,
            "-asr",
            "-quiet"
        ])

    # check the files   
    assert check_satute_files_dir(dest_dir_path, categories, alpha, asr)
    

def test_3(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test Two Without IQ-TREE File"

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
    assert excinfo.value.code == 2  # assuming exit code 1 for failure


def test_4(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test Two Without IQ-TREE File"

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
                '_bla_bla_bla',
                "-alpha",
                alpha,
                "-asr",
                "-quiet"
            ])
    # Verify the exit code if needed
    assert excinfo.value.code == 2  # assuming exit code 1 for failure    


def test_5(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test 5 MSA + DIR"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    # Satute run    
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                '-dir',
                str(dest_dir_path.absolute()),
                '-msa',
                os.path.join(dest_dir_path, "toy_example_ntaxa_7_run_5-alignment.phy"),                
                "-alpha",
                alpha,
                "-asr",
                "-quiet"
            ])
    # Verify the exit code if needed
    assert excinfo.value.code == 1  # assuming exit code 1 for failure    
    
def test_5(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test 5 MSA+DIR"

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)
    iq_tree_run_file_list = list_filenames_in_directory(source_path)
        
    # copy msa file to output directory
    copy_files_to_dest_dir(source_path, dest_dir_path, iq_tree_run_file_list)

    categories = []
    alpha = str(0.05)
    asr = True
    
    # Satute run    
    with pytest.raises(SystemExit) as excinfo:
        run_satute(
            [
                '-dir',
                str(dest_dir_path.absolute()),
                '-msa',
                os.path.join(dest_dir_path, "toy_example_ntaxa_7_run_5-alignment.phy"),                
                "-alpha",
                alpha,
                "-asr",
                "-quiet"
            ])
    # Verify the exit code if needed
    assert excinfo.value.code == 1  # assuming exit code 1 for failure    
    
def test_6(dir_path, iqtree, python, satute):
    source_path  = dir_path[0] 
    suffix = "Directory Test 6 dir option edge"

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

def test_7(dir_path_categories, iqtree, python, satute):
    source_path  = dir_path_categories[0] 
    suffix = "Directory Test 7 With Rate Categories And Edge"

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


def test_8(dir_path_categories, iqtree, python, satute):
    source_path  = dir_path_categories[0] 
    suffix = "Directory Test 8 With Rate Categories And Edge"

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
    
    
def test_9(dir_path_categories, iqtree, python, satute):
    source_path  = dir_path_categories[0] 
    suffix = "Directory Test 9 Without Siteprob file"

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
    assert excinfo.value.code == 1  # assuming exit code 1 for failure    
