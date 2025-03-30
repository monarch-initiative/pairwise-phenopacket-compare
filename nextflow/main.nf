#!/usr/bin/env nextflow

include { get_env } from './modules/get_env.nf'
include { download_and_setup } from './modules/download_and_setup.nf'
include { compute_pairwise_comps } from './modules/pairwise_comparisons.nf'
include { cluster_and_plot } from './modules/cluster_and_plot.nf'


// Phenologs calculation parameters 
params.release = "0.1.24"
params.cpu_cores = 10 // Not actuall used (any more... config takes care of this)
params.sub_sample = false
params.sim_metric = "aic"


workflow {

    Channel.value(params.release).set{ release }
    Channel.value(params.sub_sample).set{ sub_sample }
    Channel.value(params.sim_metric).set{ sim_metric }

    // Setup data environment
    get_env()
    download_and_setup(get_env.out.env_path, release)

    // Compute pairwise comparisons
    compute_pairwise_comps(get_env.out.env_path, 
                           download_and_setup.out.project_path,
                           release,
                           sub_sample)

    // Cluster and plot
    cluster_and_plot(get_env.out.env_path, 
                     download_and_setup.out.project_path,
                     compute_pairwise_comps.out.comp_sig,
                     release,
                     sub_sample,
                     sim_metric)
}