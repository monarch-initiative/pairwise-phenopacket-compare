# Imports
import os
import argparse
import scipy
import pickle
import gzip
import copy
import multiprocessing as mp
import numpy as np
import pandas as pd
from collections import Counter
from typing import List, Optional
from pydantic import BaseModel

# Data vis
import matplotlib.pyplot as plt
import xarray

# Pheval phenopacket utils, Ontology tools
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil
from semsimian import Semsimian



class Patient(BaseModel):
    sample_name: str
    phenopacket_path: str
    
    # Pulled from phenopacket
    phenotype_ids: list
    phenotype_count: int
    
    disease_name: str
    disease_id: str
    
    gene_symbol: str
    gene_id: str
    
    # For sssom mappings
    disease_mapped: Optional[str] = None
    gene_mapped: Optional[str] = None


def read_sssom_to_lookup(fpath, exact_match=True):
    """
    Assumes monarch initiative header style of # and then first line is header
    """
    
    sssom_map = {}
    with open(fpath, 'r') as infile:
        header = 0
        cols = {}
        for line in infile:

            
            line = line.strip('\r').strip('\n')
            
            # Weird thing happening with the monarch sssom files (new line character most likely as the final line)
            if len(line) == 0:
                break
                
            if line[0] == "#":
                continue
            
            header += 1
            cols = line.split('\t')
            if header == 1:
                col_inds = {v:i for i,v in enumerate(cols)}
                continue
            
            # Our actual data here
            map_key = cols[col_inds["object_id"]]
            map_val = cols[col_inds["subject_id"]]
            map_type = cols[col_inds["predicate_id"]]
            
            # Only deal with exact matches
            if map_type != "skos:exactMatch":
                continue
            
            sssom_map.update({map_key:map_val})
    
    val_count = len(set(list(sssom_map.values())))
    print("- {} unique terms mapping to {} unique terms...".format(format(len(sssom_map), ','), format(val_count, ',')))
    return sssom_map


def gather_phenopacket_paths(pheno_store_release_dir):
    """
    Assumes phenopacketstore release structure (directory per set of phenopackets)
    """
    
    phen_paths = {}
    dir_count = 0
    for fname in os.listdir(pheno_store_release_dir):
        fpath = os.path.join(pheno_store_release_dir, fname)
        if not os.path.isdir(fpath):
            continue
        
        dir_count += 1
        for pname in os.listdir(fpath):
            if not pname.endswith(".json"):
                continue
            
            sname = pname.replace(".json", '')
            phen_path = os.path.join(fpath, pname)
            if sname not in phen_paths:
                phen_paths.update({sname:phen_path})
            else:
                print("- ERROR - {} found multiple times... Current schema needs adjusting".format(sname))
    
    print("- {} Phenopackets found across {} directories...".format(format(len(phen_paths), ','), dir_count))
    return phen_paths


def phenopacket_paths_to_data(sample_path_dict):
    
    # Return data structure and function variables
    pobjs = {}
    multi_gene, multi_dis = {}, {}
    processed, total_samples = 0, len(sample_path_dict)
    
    # Loop through all phenopackets, extract / map data and return Patient objects
    for sname, spath in sample_path_dict.items():
        
        phenopacket_util = PhenopacketUtil(phenopacket_reader(spath))
        observed_phenotypes = phenopacket_util.observed_phenotypic_features()
        phenotype_ids = [observed_phenotype.type.id for observed_phenotype in observed_phenotypes]
        phenotype_count = len(phenotype_ids)
        
        # Default is to take first term 
        # Will display how many multi terms there are at end... should be few if not none)
        dis_obj = phenopacket_util.diseases()[0]
        dis_name, dis_id = dis_obj.term.label, dis_obj.term.id
        
        gene_obj = phenopacket_util.diagnosed_genes()[0]
        gene_symbol, gene_id = gene_obj.gene_symbol, gene_obj.gene_identifier
        
        if len(phenopacket_util.diagnosed_genes()) > 1:
            multi_gene.update({sname:''})
        
        if len(phenopacket_util.diseases()) > 1:
            multi_dis.update({sname:''})
        
        #print(sname)
        #print(dis_name, dis_id)
        #print(gene_symbol, gene_id)

        pobjs.update({sname:Patient.model_validate({"sample_name":sname,
                                                    "phenopacket_path":copy.copy(spath),
                                                    "phenotype_ids":phenotype_ids,
                                                    "phenotype_count":phenotype_count,
                                                    "disease_name":dis_name,
                                                    "disease_id":dis_id,
                                                    "gene_symbol":gene_symbol,
                                                    "gene_id":gene_id})})
        
        processed += 1
        if processed % 1_000 == 0:
            print("- {}/{} phenopackets read into memory...".format(format(processed, ','),
                                                                    format(total_samples, ',')))
        
    
    print("- Multi gene diagnosis phenopackets found {}...".format(len(multi_gene)))
    print("- Multi disease diagnosis phenopackets found {}...".format(len(multi_dis)))
    print("- {} Phenopackets information read into memory...".format(format(len(pobjs), ',')))
    return pobjs


def filter_non_zero_data(phen_data, sub_sample=False):
    # Remove zero phenotype count samples
    removed = 0
    for k in list(phen_data.keys()):
        if phen_data[k].phenotype_count == 0:
            del phen_data[k]
            removed += 1       
    print("- {} samples removed with zero phenotypes...".format(removed))
    
    # Default is no subsampleing
    if sub_sample != False:
    
        # Subsample (to make full pipeline connection easier)
        sub_samp = 0
        removed = 0
        for k in list(phen_data.keys()):
            if sub_samp >= sub_sample:
                del phen_data[k]
                removed += 1

            sub_samp += 1
        print("- {} samples removed with zero phenotypes...".format(removed))
        print("- {} samples remaining...".format(len(phen_data)))
    
    return phen_data
    
    
def generate_pairwise_comparison_args(sample_data):
    
    keys = list(sample_data.keys())
    comps = []
    ind = 0
    for i,k in enumerate(keys[0:-1]):
        p_ids = sample_data[k].phenotype_ids
        for k2 in keys[i+1:]:
            p_ids2 = sample_data[k2].phenotype_ids
            comps.append([k, k2, p_ids, p_ids2, ind])
            ind += 1
    
    print("- {} Pairwise comparison terms generated...".format(format(len(comps), ',')))
    print("- {} Labels generated...".format(format(len(keys), ',')))
    return comps, keys


def semsim_termset_compare(hp_db_path, input_args):
    
    # Load in data
    semsim_obj = Semsimian(spo=None, resource_path=hp_db_path)

    # Create output datastructure
    sim_metrics = ["jaccard_similarity", "ancestor_information_content", "phenodigm_score"]
    output_data = {"comp_name":[],
                   "comp_index":[],
                   "jaccard_similarity":[],
                   "ancestor_information_content":[],
                   "phenodigm_score":[]}
    
    # For each set of input arguments, compute jaccar, ic, and phenodigm sim scores and add to our output dict
    for i, argset in enumerate(input_args):
        
        termset1, termset2 = argset[2], argset[3]
        sname1, sname2 = argset[0], argset[1]
        comp_index = argset[4]
        output_data["comp_name"].append(tuple([sname1, sname2]))
        output_data["comp_index"].append(comp_index)
        
        for sim in sim_metrics:
            score = semsim_obj.termset_comparison(set(termset1), set(termset2), score_metric=sim)
            output_data[sim].append(score)

        if i % 1000 == 0:
            print("- processed {}/{}".format(i+1, format(len(input_args), ',')))

    return output_data


def merge_and_sort_parallel_output(output, labels, outpath=False):
    
    # Merge results from previous step into single data structure
    # Each element in the output list is a dictionary {sival:[pval,pval,...], }
    merged_data = {}
    for p in output:
        for k,v in p.items():
            if k not in merged_data:
                merged_data.update({k:[]})

            if type(v) != type([]):
                v = list(v)
            merged_data[k] += v
    
    for k,v in merged_data.items():
        print(k, len(v))
    # Sort all data based on the "comp_index"... 
    # This allows us to recapitulate the input "collapsed" pairwise comparison matrix
    for mkey in list(merged_data.keys()):
        print(mkey, "merging...")

        if mkey == "comp_index":
            continue
        
        # Must avoid combining strings with integers for array notation
        if mkey == "comp_name":
            cmat = sorted([[str(v),str(ind)] for v,ind in zip(merged_data[mkey], merged_data["comp_index"])], key=lambda x: int(x[1]))
        
        else:
            cmat = sorted([[v,ind] for v,ind in zip(merged_data[mkey], merged_data["comp_index"])], key=lambda x: x[1])
        
        cmat = np.asarray(cmat).T[0]
        merged_data[mkey] = cmat
    
    merged_data["comp_index"] = sorted(merged_data["comp_index"])
    merged_data.update({"labels":labels})
    if outpath != False:
        pickle.dump(merged_data, gzip.open(outpath, "wb"))
    
    return merged_data



def normalize_condensed_matrix(condensed_matrix, to_dist=True, max_val=0):
    """
    This function is a general function to convert back and forth from a distance to a similarity matrix
       Input = adjacencyMatrix and corresponding Bin object list (we use the preComputed rowSum attribute of the binObjects to speed up the process )
       Output = transformed adjacencyMatrix
    """
    
    if type(condensed_matrix) == type([]):
        condensed_matrix = np.asarray(condensed_matrix)
    
    norm_val = max(condensed_matrix)
    if to_dist == True:
        condensed_matrix = condensed_matrix / norm_val
    else:
        condensed_matrix = condensed_matrix * max_val
    
    return condensed_matrix, norm_val


def reorderMatrix(matrix,binList,newOrder):
    """
    function to reorder the rows and columns of our matrix according to new indices (also the list of Bin objects)
       Input = adjacencyMatrix, list of Bin objects, indices of newOrder
       Output = reordered matrix and list of Bin objects
    """ 
    matrix = matrix[:, newOrder][newOrder]
    binList = [ binList[i] for i in newOrder ]
    return matrix,binList


def averageClusterNodes(condensedMatrix,nodeLabels,noPlot=True):
    '''Performs UPGMA (average) heirarchical clustering for a distance transformed adjacency matrix
       The actual matrix passed in is not altered in any way. It must be roered in another function
       Input = distance transformed adjacencyMatrix
               labels of the nodes of the dendrogram that is produced
       Output = dendrogram object produced by scipy.cluster.dendrogram'''
    
    clustMan = scipy.cluster.hierarchy.average(condensedMatrix)
    print("- Data formatted for clustering...")
    
    if noPlot == False:
        fig = plt.figure()
        fig.set_size_inches(32,8)
        dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='ascending')
        plt.show()
    else:
        dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=noPlot,get_leaves=True,count_sort ='ascending')
    
    print("- Heirarchical clustering complete...")
    return dendro


def dendrogramLeafOrder_toFile(dendrogramObj,outFile):
    '''This is to simply write out the order of the nodes of the dendrogram produced by the averageClusterNodes function
       Input = dendrogram object from scipy, full file path to output file'''
    outFile = open(outFile,'w')
    for i,l in enumerate(dendrogramObj['ivl']):
        if i != (len(dendrogramObj['ivl'])-1):
            outFile.write(l+'\t'+str(dendrogramObj['leaves'][i])+'\n')
        else:
            outFile.write(l+'\t'+str(dendrogramObj['leaves'][i]))
    ###############
    outFile.close()

    
def readDengrogramLeavesFromFile(dendrogramFile):
    '''This simply reads in the output file from the dendrogramLeafOrder_toFile
       Input = full file path to dendrogramFile from dendrogramLeafOrder_toFile function
       Output = a dictionary with the same notation as the scipy dendrogram object, but only "ivl" and "leaves" are used" '''
    inFile = open(dendrogramFile,'r')
    dendoLeaves = { 'ivl':[],'leaves':[] }
    for line in inFile:
        cols = line.strip('\r').strip('\n').split('\t')
        dendoLeaves['ivl'].append(cols[0])
        dendoLeaves['leaves'].append(int(cols[-1]))
    ##############
    inFile.close()
    return dendoLeaves


def plotContactMap(adjMat,
                   wInches=32,
                   hInches=32,
                   lP=1,
                   hP=98,
                   reverseColorMap='_r',
                   showPlot=False,
                   savePlot=False,
                   title=False,
                   titleSuffix=False):
    
    if type(adjMat) != type(np.asarray([])):
        print("- Converting input type {} --> {} for plotting functionality".format(type(adjMat), type(np.asarray([]))))
        adjMat = np.asarray(adjMat)

    print("- Attempting to plot array / matrix of shape {}".format(adjMat.shape))
    matrixLength = len(adjMat)
    fig,ax = plt.subplots()
    fig.set_size_inches(wInches,hInches)
    xarray.plot.pcolormesh(xarray.DataArray(adjMat[::-1]),
                         ax=ax,
                         add_colorbar=False,
                         cmap="plasma"+reverseColorMap, #plasma_r for the reverse which is best for distance matricies. Normal or '' is best for similarity matrices
                          
                         robust=True,
                         vmin=np.percentile(adjMat,lP),
                         vmax=np.percentile(adjMat,hP))


    #tickDist = len(adjMat) / tickCount
    #xTicks,tickMan = [0],0
    #for i in range(0,tickCount-1):
    #    tickMan += tickDist
    #    xTicks.append(tickMan)

    #xTicks.append(len(adjMat))
    #ax.set_xticks(xTicks)
    #ax.set_xticklabels([ str(int((t*resolution)/1000000))+" Mb" for t in xTicks ],size=18)
    #ax.set_xlabel('')
    #xTicks.pop(0)
    #ax.set_yticks(xTicks)
    #ax.set_yticklabels([ str(int((t*resolution)/1000000))+" Mb" for t in xTicks ],size=18)
    #ax.set_ylabel('')

    if title != False:
        if titleSuffix != False:
            title = title+titleSuffix
        ax.set_title(title,size=25)

    if savePlot != False:
        plt.savefig(savePlot)

    if showPlot != False:
        plt.show()
    plt.close(fig)


def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Phenopacket Pairwise Similarity...')
        parser.add_argument("-p","--project_dir", help="Top most project directory", required=True, type=str, default=None)
        parser.add_argument("-r","--release", help="Phenopackets store release version (0.1.24, 0.1.23, etc...)", required=True, type=str, default=None)
        parser.add_argument("-c","--num_proc", help="number of cpu cores to use", required=False, type=int, default=2)
        parser.add_argument("-s","--subsample", help="subsample data for testing", required=False, default=False)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    #all_phens_dir = "/Users/ao33/Desktop/PHENOPACKETS_STORE/0.1.24/"
    all_phens_dir = os.path.join(args.project_dir, "phenopackets_store", args.release)
    mondo_sssom_path = os.path.join(args.project_dir, "ontology_files", "mondo.sssom.tsv")
    gene_sssom_path = os.path.join(args.project_dir, "ontology_files", "gene_mappings.sssom.tsv")
    hp_db_path = os.path.join(args.project_dir, "ontology_files", "hp.db")


    if args.subsample == False:
        outname = "{}_results.pkl.gz".format(args.release)
    else:
        outname = "{}_results_subsample_{}.pkl.gz".format(args.release, args.subsample)

    res_path = os.path.join(args.project_dir, "results", outname)

    # Create mappings
    mondo_map = read_sssom_to_lookup(mondo_sssom_path)
    gene_map = read_sssom_to_lookup(gene_sssom_path)

    # Load phenopacket data to memory
    phen_paths = gather_phenopacket_paths(all_phens_dir)
    phen_data = phenopacket_paths_to_data(phen_paths)

    # Filter non-zero phenotype count samples
    phen_data = filter_non_zero_data(phen_data, sub_sample=args.subsample)
    comp_args, labs = generate_pairwise_comparison_args(phen_data)

    # Divide workload for parallel processing
    div_args = divide_workload(comp_args, num_proc=args.num_proc)
    div_semsim_objs = [copy.copy(hp_db_path) for i in range(0, args.num_proc)]

    # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
    output = mp.Queue()
    pool = mp.Pool(processes=args.num_proc)
    results = [pool.apply_async(semsim_termset_compare, args=(sem, inargs,)) for sem, inargs in zip(div_semsim_objs, div_args)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    print("- Semantic similarity comparisons completed...")
    print("- Combining and writing results to {}...".format(res_path))
    comp_data = merge_and_sort_parallel_output(output=output, labels=labs, outpath=res_path)