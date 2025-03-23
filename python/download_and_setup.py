# General imports
import os
import argparse
import requests
import tarfile
import pickle
import zipfile
import gzip
import shutil
 

def download_file_url(url: str, outdir: str, extract_gz: bool = False, overwrite: bool = False, extract_zip: bool = False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/"
    """
    
    # Download and write file
    filename = os.path.join(outdir, url.split("/")[-1])
    fname = filename.split("/")[-1]

    if overwrite == True:
        if os.path.isfile(filename):
            print("- Warning, file {} already exists... Set overwrite to True to download and replace")
            return

    with open(filename, "wb") as f:
        r = requests.get(url)
        f.write(r.content)
    
    # Extract gzip
    if extract_gz != False:
        outpath = os.path.join(outdir, fname.replace(".gz", ""))
        with gzip.open(filename, 'rb') as f_in:
            with open(outpath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    elif extract_zip != False:
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(outdir)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Performs phenologs setup steps. Makes relevant directory structures, \
                                                      downloads monarch kg, and upacks the kg files. Will not overwrite \
                                                      existing data or monarch-kg data if already exists within specified \
                                                      projece_dir')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        parser.add_argument("-r", "--release", help="Phenopackets store release version (0.1.24, 0.1.23, etc...)", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # Ontology / mappings download paths, and phenopackets store
    onto_dir = os.path.join(args.project_dir, "ontology_files")
    pheno_dir = os.path.join(args.project_dir, "phenopackets_store")

    mondo_path = os.path.join(onto_dir, "mondo.sssom.tsv")
    gene_path = os.path.join(onto_dir, "gene_mappings.sssom.tsv")
    hp_path = os.path.join(onto_dir, "hp.db")


    # Project top level data directories / structure 
    project_dirs = ["phenopackets_store",
                    "ontology_files",
                    "results"]

    # Create base project directory
    if not os.path.isdir(args.project_dir):
        print("- Creating project directory(ies) at {}".format(args.project_dir))
        os.makedirs(args.project_dir, exist_ok=True)

    # Create dataset directories
    for pdir in project_dirs:
        os.makedirs(os.path.join(args.project_dir, pdir), exist_ok=True)

    # Download and upack hp.db file
    if not os.path.isfile(hp_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking hp.db to {}".format(hp_path))
        URL = 'https://s3.amazonaws.com/bbop-sqlite/hp.db.gz'
        download_file_url(URL, onto_dir, extract_gz=True, extract_zip=False, overwrite=False)
        print("- Download and upacking of hp succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(hp_path))
    

    # Download and upack gene sssom file
    if not os.path.isfile(gene_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking gene sssom to {}".format(gene_path))
        URL = 'https://data.monarchinitiative.org/mappings/latest/gene_mappings.sssom.tsv'
        download_file_url(URL, onto_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of genes succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(gene_path))
    

    # Download and upack mondo sssom file
    if not os.path.isfile(mondo_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking mondo sssom to {}".format(mondo_path))
        URL = 'https://data.monarchinitiative.org/mappings/latest/mondo.sssom.tsv'
        download_file_url(URL, onto_dir, extract_gz=False, extract_zip=False, overwrite=False)
        print("- Download and upacking of mondo succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(mondo_path))
    

    # Download phenopacket store data
    pheno_release_dir = os.path.join(pheno_dir, args.release)
    if not os.path.isdir(pheno_release_dir):
        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking mondo sssom to {}".format(mondo_path))
        URL = 'https://github.com/monarch-initiative/phenopacket-store/releases/download/{}/all_phenopackets.zip'.format(args.release)
        download_file_url(URL, pheno_dir, extract_gz=False, extract_zip=True, overwrite=False)
        print("- Download and upacking of phenopackets succesfull...")
    else:
        print("- Skipping hp download... File already exists at {}".format(pheno_release_dir))
    
    
    print("- Setup complete...")