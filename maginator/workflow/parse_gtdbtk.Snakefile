import os

WD = config['wd']
PARAMS = config['params']
VAMB = config['vamb']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes95_cluster.tsv'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers_bins_geneID.tsv')

rule parse_gtdbtk:
    input:
        os.path.join(WD, 'gtdbtk'),
        WD
    output:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes.fna'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_bac_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_ar_markers.tab')
    params:
        param_dict['annotation_prevalence'],
        VAMB
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores=1,
        memory=32,
        runtime='10:00:00'
    script:
        "../parse_gtdbtk.py"


# Get representative genes from all genes.
rule repres_genes:
    input:
        os.path.join(WD, 'genes', 'all_genes.fna'),
    output:
        fasta = os.path.join(WD, 'genes', 'all_genes95_rep_seq.fasta'),
        tsv = os.path.join(WD, 'genes', 'all_genes95_cluster.tsv')
    resources:
        cores = 14,
        memory = 50,
        runtime = '2:00:00:00' 
    params:
        tmp_dir = os.path.join(WD, 'tmp'),
        out_prefix = os.path.join(WD, 'genes', 'all_genes95')
    conda:
        "envs/filter_gtdbtk.yaml"
    shell:
        "mmseqs easy-linclust --min-seq-id 0.95 -c 0.95 --threads {threads} {input} {params.out_prefix} {params.tmp_dir}; rm -r {params.tmp_dir};"


# Add gene clusters to GTDB-tk data
rule join:
    input:
        gtdb = os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        cluster = os.path.join(WD, 'genes', 'all_genes95_cluster.tsv')
    output:
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers_bins_geneID.tsv')
    resources:
        cores = 1,
        memory = 20,
        runtime = '24:00:00'
    shell:
        "join -1 1 -2 2 <(sort -k1,1 {input.gtdb}) <(sort -k2,2 {input.cluster}) > {output}"

