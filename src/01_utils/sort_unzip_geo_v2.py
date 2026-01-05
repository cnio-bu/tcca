import glob
import gzip
import os
import re
import shutil
import sys


def unzip_and_save(gz_file:str, where_to_save:str):
    
    with gzip.open(gz_file, "rb") as f_in:
        with open(where_to_save, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

def main(dir: str):
    
    ## populate a list of sample names based on GEO Syntax
    all_samples_files = glob.glob(os.path.join(dir + "/GSM*.gz"))
    
    ## build a sample set based on matching
    sample_name = re.compile("GSM[0-9]+")
    
    samples = [sample_name.search(f) for f in all_samples_files]
    samples = set([s.group(0) for s in samples])
    
    for sample in samples:
        
        ## check if a subdir. exists 
        dir_exists = os.path.exists(os.path.join(dir + sample))
        
        if not dir_exists:
            
            os.makedirs(os.path.join(f"{dir}/{sample}"))
            
            ## now unzip and rename
            mat_file = glob.glob(f"{dir}/{sample}_*_matrix.mtx.gz")[0]
            genes_file = glob.glob(f"{dir}/{sample}_*_features.tsv.gz")[0]
            barcodes_file = glob.glob(f"{dir}/{sample}_*_barcodes.tsv.gz")[0]
            
            unzip_and_save(mat_file, f"{dir}/{sample}/matrix.mtx")
            unzip_and_save(genes_file, f"{dir}/{sample}/genes.tsv")
            unzip_and_save(barcodes_file, f"{dir}/{sample}/barcodes.tsv")

            
    

if __name__ == "__main__":
    main(sys.argv[1])
    
    
    
    
