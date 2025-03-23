
process get_env {
    tag "get_env"
    publishDir "./", mode: 'copy'

    output:
    path "pairwise-phenopacket-compare", emit: env_path

    script:
    """
    git clone monarch_phenologs_optimized_fdr_calcs git@github.com:monarch-initiative/pairwise-phenopacket-compare.git
    cd pairwise-phenopacket-compare
    poetry config virtualenvs.in-project true
    poetry install
    """
}