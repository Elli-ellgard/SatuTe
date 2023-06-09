import os

def remove_all_except(directory, ignore_list):
    for file in os.listdir(directory):
        if file not in ignore_list:
            os.remove(os.path.join(directory, file))

import shutil

def delete_folder(folder_path):
    """Deletes the specified folder."""
    try:
        shutil.rmtree(folder_path)
        print(f"Successfully deleted the folder: {folder_path}")
    except OSError as e:
        print(f"Error: {folder_path} : {e.strerror}")


if __name__ == "__main__":
    directory = "./test/octo-kraken-msa-test/"
    keep_files = ["example.phy", "example.phy.treefile", "example.phy.tree"]
    remove_all_except(directory, keep_files)
    # delete_folder("./test/octo-kraken-msa-test/clades")
