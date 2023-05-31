import os

def remove_all_except(directory, ignore_list):
    for file in os.listdir(directory):
        if file not in ignore_list:
            os.remove(os.path.join(directory, file))

if __name__ == "__main__":
    directory = "./test/octo-kraken-msa-test/"
    keep_files = ["example.phy", "example.phy.treefile"]
    remove_all_except(directory, keep_files)

