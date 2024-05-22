from tests.scripts.satute_test_utils import (
    run_external_command,
    print_test_name, 
    create_destination_dir,
    copy_files_to_dest_dir,
    check_iqtree_files_exist,
    check_satute_files,
    print_colored_message
    
)

from tests.scripts.fixtures import *
import pytest

@pytest.mark.skip(reason='Nothing is being tested, test incomplete')
def test_siteprob(data_dir_path,iqtree, python, satute):
    source_path, msa, treefile = data_dir_path
    suffix = "SITEPROB TEST 1"
    print_test_name(suffix)

    # create output directory
    dest_dir_path = create_destination_dir(source_path, suffix)

    # copy msa file to output directory
    files_to_copy = [msa]
    copy_files_to_dest_dir(source_path, dest_dir_path, files_to_copy)

    categories = []
    alpha = str(0.05)
    asr = True
    
#     # Satute run
#     run_external_command(
#         [
#             python,
#             satute,
#             "-iqtree",
#             iqtree,
#             "-msa",
#             os.path.join(dest_dir_path, msa),
#             "-alpha",
#             alpha,
#             "-asr",
#         ]
#    )

#     # check the files   
#     if check_iqtree_files_exist(msa, dest_dir_path, ["m"] )  and check_satute_files(msa, dest_dir_path, categories, alpha, asr):
#         print_colored_message(f"{suffix} was successful", "32" )
#     else: 
#         print_colored_message(f"{suffix} failed", "31" )




