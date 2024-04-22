WD = config['wd']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(F_DIR, '{cluster}')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

SAMPLES = set()
sample_dict = {}

with open(param_dict['reads']) as f:
    for row in f:
        row = row.rstrip().split(',')
        SAMPLES.add(row[0])
        sample_dict[row[0]] = [row[1], row[2]]

rule all:
    input:
        expand(os.path.join(WD,'unknown_bugs','bams','{sample}.filtered.bam'), sample=SAMPLES)

rule concat:
    input:
        os.path.join(WD,'clusters')
    output:
        os.path.join(WD, 'unknown_bugs', 'all_mags.fasta')
    resources:
        cores = 1,
        mem_gb = 50,
        runtime = 43200 #12h in s
    shell:
        'cat {input}/*/*.fa > {output}'


rule bwa_index:
    input: 
        fasta=os.path.join(WD, 'unknown_bugs', 'all_mags.fasta')
    output: 
        index=os.path.join(WD, 'unknown_bugs', 'all_mags')
    conda:
        "envs/filter_gtdbtk.yaml" 
    resources:
        cores = 40,
        mem_gb = 188, 
        runtime = 86400 #1d in s
    shell:
        "bwa-mem2 index -p {output.index} {input.fasta}; samtools faidx {input.fasta}; touch {output.index}"


# Readmapping
rule bwa_readmap:
    input:
        index = os.path.join(WD, 'unknown_bugs', 'all_mags'),
        fasta = os.path.join(WD, 'unknown_bugs', 'all_mags.fasta'),
        fastq1 = lambda wildcards: sample_dict[wildcards.sample][0],
        fastq2 = lambda wildcards: sample_dict[wildcards.sample][1]
    output:
        bam = os.path.join(WD, 'unknown_bugs', 'bams', '{sample}.bam')
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores = 40,
        mem_gb = 188,
        runtime = 86400 #1d in s
    shell:
        """
        bwa-mem2 mem -t {resources.cores} {input.index} {input.fastq1} {input.fastq2} | \
        samtools view -T {input.fasta} -b --threads {resources.cores} > {output}
        """

rule filter_bamfile:
    input:
        os.path.join(WD,'unknown_bugs', 'bams','{sample}.bam'),
    output:
        os.path.join(WD,'unknown_bugs','bams','{sample}.filtered.bam')
    conda:
        "envs/signature_reads.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    params: 
        min_map = param_dict['min_map'],
        min_iden = param_dict['min_identity'],
        min_len = param_dict['min_length'],
    shell:
        """
        msamtools filter -b -l {params.min_len} -p {params.min_iden} -z {params.min_map} --besthit {input} > {output}
        """  
