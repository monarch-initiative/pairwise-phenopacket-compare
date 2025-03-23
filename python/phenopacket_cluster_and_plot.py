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


def averageClusterNodes(condensedMatrix,nodeLabels, showPlot=True, savePlot=False):
    '''Performs UPGMA (average) heirarchical clustering for a distance transformed adjacency matrix
       The actual matrix passed in is not altered in any way. It must be roered in another function
       Input = distance transformed adjacencyMatrix
               labels of the nodes of the dendrogram that is produced
       Output = dendrogram object produced by scipy.cluster.dendrogram'''
    
    clustMan = scipy.cluster.hierarchy.average(condensedMatrix)
    print("- Data formatted for clustering...")

    fig = plt.figure()
    fig.set_size_inches(32,8)
    dendro = scipy.cluster.hierarchy.dendrogram(clustMan,labels=nodeLabels,leaf_rotation=90,no_plot=showPlot,get_leaves=True,count_sort ='ascending')
    
    if savePlot != False:
        plt.savefig(savePlot)
    
    if showPlot != False:
        plt.show()
    
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






if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Phenopacket Pairwise Similarity...')
        parser.add_argument("-p","--project_dir", help="Top most project directory", required=True, type=str, default=None)
        parser.add_argument("-r","--release", help="Phenopackets store release version (0.1.24, 0.1.23, etc...)", required=True, type=str, default=None)
        parser.add_argument("-s","--subsample", help="subsample data for testing", required=False, type=int, default=False)
        parser.add_argument("-m","--sim_metric", help="similarity metric to use", required=False, type=str, default="jaccard_similarity")
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Grab relevant filename and load data
    if (args.subsample == False) or (args.subsample == "False") or (args.subsample == "false"):
        outname = "{}_results.pkl.gz".format(args.release)
    else:
        args.subsample = int(args.subsample)
        outname = "{}_results_subsample_{}.pkl.gz".format(args.release, args.subsample)

    inpath = os.path.join(args.project_dir, "results", outname)
    comp_data = pickle.load(gzip.open(inpath, "rb"))
    labs = comp_data["labels"]
    print("- Loaded comparison data from {}. Clustering data...".format(inpath))

    # Output figure paths
    dendro_path = os.path.join(args.project_dir, 
                               "results", 
                               outname.replace(".pkl.gz", "_{}_upgma_dendro.png".format(args.sim_metric)))

    heatmap_path = os.path.join(args.project_dir, 
                                "results", 
                                outname.replace(".pkl.gz", "_{}_upgma_heatmap.png".format(args.sim_metric)))
    
    # Similarity metric mapping
    sim_metrics = {"ps": "phenodigm_score",
                   "js": "jaccard_similarity",
                   "aic": "ancestor_information_content"}
    sim_metric = sim_metrics[args.sim_metric]


    # Cluster calculations
    dendro = averageClusterNodes(normalize_condensed_matrix(comp_data[sim_metric])[0], labs, showPlot=False, savePlot=dendro_path)
    expanded_matrix = scipy.spatial.distance.squareform(comp_data[sim_metric], checks=False)
    expanded_matrix, labs = reorderMatrix(expanded_matrix, labs,dendro["leaves"])
    collapsed_matrix = scipy.spatial.distance.squareform(expanded_matrix, checks=False)
    print("- Data clustered. Plotting heatmap and saving results to {}...".format(heatmap_path))

    plotContactMap(expanded_matrix,
                   wInches=32,
                   hInches=32,
                   lP=1,
                   hP=98,
                   reverseColorMap='_r',
                   showPlot=False,
                   savePlot=heatmap_path,
                   title=False,
                   titleSuffix=False)