#!/usr/bin/env nextflow

process compute_pairwise_comps {
    tag 'compute_pairwise_comps'
    publishDir "./", mode: 'copy'

    input:
    path env_dir
    path data_dir
    val release
    val sub_sample

    output:
    path data_dir, emit: project_path
    val "done", emit: comp_sig

    script:
    """
    source ${env_dir}/.venv/bin/activate
    python pairwise-phenopacket-compare/python/phenopacket_pairwise_similarity.py -p ./pairwise-comparisons-data \
                                                                                  -r ${release} \
                                                                                  -c ${task.cpus} \
                                                                                  -s ${sub_sample} \
    """
}