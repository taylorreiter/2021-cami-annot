import pandas as pd

TMPDIR = "/scratch/tereiter"

m = pd.read_csv("inputs/CAMI_low_genomes_and_contigs.tsv", sep = "\t", header = 0)
SOURCE_GENOMES = m['BINID'].unique().tolist()
CONTIGS = m['contig'].unique().tolist()
SOURCE_GENOME_AND_CONTIG= m['genome_and_contig'].unique().tolist()

DATASETS = ['CAMI_low']

rule all:
    input:
        expand("outputs/bowtie2/{dataset}_unmapped_R2.fq", dataset = DATASETS),
        ancient(expand("outputs/eggnog_source_genomes/{source_genome}.emapper.annotations", source_genome = SOURCE_GENOMES)),
        "outputs/outputs/gs_read_annotations/CAMI_low_gs_read_annotations.tsv"

#############################################################
## Obtaining data
#############################################################

rule download_CAMI:
    output: "inputs/CAMI_low.tar"
    threads: 1
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    shell:'''
    wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100344/ChallengeDatasets.dir/CAMI_low.tar
    '''

rule decompress_CAMI:
    input: "inputs/CAMI_low.tar"
    threads: 1
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    output: 
        "inputs/CAMI_low/RL_S001__insert_270.fq.gz",
        "inputs/CAMI_low/source_genomes_low.tar.gz",
        "inputs/CAMI_low/gsa_mapping.binning",
        "inputs/CAMI_low/gs_read_mapping.binning.gz",
    shell:'''
    tar xvf {input} -C inputs/
    '''

#################################################################
## Identifying unassembled reads
#################################################################

rule fastp:
    """
    The CAMI I and CAMI II challenge both showed that "using read quality trimming or error correction software, such as ...Fastp... impoved assembly quality." https://doi.org/10.1101/2021.07.12.451567
    Set minimum read length to megahit default minimum kmer length, k = 21
    """
    input: "inputs/{dataset}/RL_S001__insert_270.fq.gz"
    output: 
        r1 = "outputs/fastp/{dataset}_R1.fastp.fq.gz",
        r2 = "outputs/fastp/{dataset}_R2.fastp.fq.gz",
        json = "outputs/fastp/{dataset}.json",
        html = "outputs/fastp/{dataset}.html"
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 8
    benchmark: "benchmarks/fastp_{dataset}.txt"
    conda: "envs/fastp.yml"
    shell:'''
    fastp -i {input} --interleaved_in -o {output.r1} -O {output.r2} -q 4 -j {output.json} -h {output.html} -R {wildcards.dataset} -l 21 -c -w {threads}
    '''

rule assemble:
    """
    The CAMI II challenge indicated that megahit and metaSPAdes performed approximately equally for assembly accuracy and strain recall. Given that megahit requires less ram/runtime, start with that assembler.
    """
    input:
        r1 = "outputs/fastp/{dataset}_R1.fastp.fq.gz",
        r2 = "outputs/fastp/{dataset}_R2.fastp.fq.gz",
    output: "outputs/megahit/{dataset}.contigs.fa"
    params: outdir = lambda wildcards: "outputs/megahit/" + wildcards.dataset + "_tmp/"
    resources:
        mem_mb = 32000,
        tmpdir=TMPDIR
    threads: 1
    benchmark: "benchmarks/megahit_{dataset}.txt"
    conda: "envs/megahit.yml"
    shell:'''
    megahit -1 {input.r1} -2 {input.r2} -o {params.outdir} --out-prefix {wildcards.dataset} --min-contig-len 500
    mv {params.outdir}/{wildcards.dataset}.contigs.fa {output}
    '''

rule index_assembly:
    """
    The CAMI I challenge used bowtie2 with --end-to-end parameter to assess the number of reads that mapped back to the assembly.
    """
    input: "outputs/megahit/{dataset}.contigs.fa"
    output: "outputs/bowtie2_index/{dataset}.1.bt2"
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    params: prefix = lambda wildcards: "outputs/bowtie2_index/" + wildcards.dataset
    benchmark: "benchmarks/bowtie2_index_{dataset}.txt"
    conda: "envs/bowtie2.yml"
    threads: 1
    shell:'''
    bowtie2-build {input} {params.prefix}
    '''

rule map_reads_to_assembly:
    """
    The CAMI I challenge used bowtie2 with --end-to-end parameter to assess the number of reads that mapped back to the assembly.
    """
    input:
        index="outputs/bowtie2_index/{dataset}.1.bt2",
        r1 = "outputs/fastp/{dataset}_R1.fastp.fq.gz",
        r2 = "outputs/fastp/{dataset}_R2.fastp.fq.gz",
    output: "outputs/bowtie2/{dataset}.bam"
    params: prefix = lambda wildcards: "outputs/bowtie2_index/" + wildcards.dataset
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/bowtie2_{dataset}.txt"
    conda: "envs/bowtie2.yml"
    threads: 8
    shell:'''
    bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} -p {threads} --end-to-end | \
    samtools view -Sbh --threads {threads} - > {output}
    '''

rule identify_unmapped_reads:
    input: "outputs/bowtie2/{dataset}.bam"
    output: "outputs/bowtie2/{dataset}_unmapped.bam"
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/samtools_f4_{dataset}.txt"
    threads: 1
    conda: 'envs/bowtie2.yml'
    shell:'''
    samtools view -b -f 4 {input} > {output}
    '''

rule convert_unmapped_reads_to_fastq:
    input:"outputs/bowtie2/{dataset}_unmapped.bam"
    output: 
        r1="outputs/bowtie2/{dataset}_unmapped_R1.fq",
        r2="outputs/bowtie2/{dataset}_unmapped_R2.fq",
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/samtools_fastq_unmapped_{dataset}.txt"
    conda: 'envs/bowtie2.yml'
    threads: 1
    shell:'''
    samtools fastq -N -1 {output.r1} -2 {output.r2} {input}
    '''

############################################################
## Generating silver-standard annotations of source genomes
############################################################

rule decompress_source_genomes:
    input: "inputs/CAMI_low/source_genomes_low.tar.gz"
    output: expand("inputs/CAMI_low/source_genomes/{source_genome}.fna", source_genome = SOURCE_GENOMES)
    threads: 1
    resources: 
        mem_mb = 500,
        tmpdir=TMPDIR
    shell:'''
    tar xvf {input} -C inputs/CAMI_low
    # these next few lines are a very poor idea, but the source genomes
    # 1) have two separate file endings (fasta, fna);
    # 2) the file prefix does not match the genome identifier used by the rest of the
    #    documents in CAMI low.
    # The following lines 
    # 1) give all source genomes the same file ending (.fna)
    # 2) unify the source genome file prefixes with those used in the rest of the 
    #    CAMI low documents.  
    for infile in inputs/CAMI_low/source_genomes/*fasta 
    do
    bn=$(basename $infile .fasta)
    mv $infile ${{bn}}.fna
    done
    
    for infile in *.gt1kb.fna
    do
    bn=$(basename $infile .gt1kb.fna)
    mv $infile ${{bn}}.fna
    done

    mv 1139_AG_run158_run162.final.scaffolds.fna 1139_AG.fna
    mv 1220_AD_run172_run176.final.scaffolds.fna 1220_AD.fna
    mv 1220_AJ_run172_run176_run188.final.scaffolds.fna 1220_AJ.fna
    mv 1285_BH_run189.final.scaffolds.fna 1285_BH.fna
    mv 1286_AP_run191_run197.final.scaffolds.fna 1286_AP.fna
    mv 1365_A_run201.final.scaffolds.fna 1365_A.fna
    '''

# uncomment to download db; use db that was downloaded for another project
# rule download_bakta_db:
#    output: "inputs/bakta_db/db/version.json"
#    threads: 1
#    resources: mem_mb = 4000
#    conda: "envs/bakta.yml"
#    shell:'''
#    bakta_db download --output {output}
#    '''

rule bakta_source_genomes:
    input: 
        fna=ancient("inputs/CAMI_low/source_genomes/{source_genome}.fna"),
        db="/home/tereiter/github/2021-orpheum-refseq/inputs/bakta_db/db/version.json"
    output: 
        "outputs/bakta_source_genomes/{source_genome}.faa",
        "outputs/bakta_source_genomes/{source_genome}.gff3",
        "outputs/bakta_source_genomes/{source_genome}.fna",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    benchmark: "benchmarks/bakta_{source_genome}.txt"
    conda: 'envs/bakta.yml'
    params: 
        dbdir="~/github/2021-orpheum-refseq/inputs/bakta_db/db/",
        outdir = 'outputs/bakta_source_genomes/',
    threads: 1
    shell:'''
    bakta --db {params.dbdir} --prefix {wildcards.source_genome} --output {params.outdir} \
        --locus-tag {wildcards.source_genome} --keep-contig-headers {input.fna}
    '''

rule remove_fasta_from_gff:
    input: gff="outputs/bakta_source_genomes/{source_genome}.gff3",
    output: no_fasta="outputs/bakta_source_genomes/{source_genome}_no_fasta.gff"
    resources: 
        mem_mb = 2000 ,
        tmpdir= TMPDIR
    threads: 1
    shell:'''
    sed '/^##FASTA$/,$d' {input} > {output}
    '''

rule download_eggnog_db:
    output: "inputs/eggnog_db/eggnog.db"
    threads: 1   
    resources: 
        mem_mb = 4000,
        tmpdir=TMPDIR
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir inputs/eggnog_db
    '''

rule eggnog_annotate_source_genomes:
    input: 
        faa = "outputs/bakta_source_genomes/{source_genome}.faa",
        db = 'inputs/eggnog_db/eggnog.db'
    output: "outputs/eggnog_source_genomes/{source_genome}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        tmpdir=TMPDIR
    threads: 8
    params: 
        outdir = "outputs/eggnog_source_genomes/",
        dbdir = "inputs/eggnog_db"
    shell:'''
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.source_genome} \
       --output_dir {params.outdir} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
    '''

#########################################################################
## Determine read base pair coordinates in source_genomes
#########################################################################

rule write_read_names_by_genome_and_contig:
    input: gs = "inputs/CAMI_low/gs_read_mapping.binning.gz"
    output: expand("outputs/gs_read_mapping/{source_genome_and_contig}.txt", source_genome_and_contig = SOURCE_GENOME_AND_CONTIG) 
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 128000,
        tmpdir=TMPDIR
    threads: 1
    script: "scripts/write_read_names_by_genome_and_contig.R"

rule split_reads_by_genome_and_contig:
    input: 
        lst="outputs/gs_read_mapping/{source_genome}-{contig}.txt",
        fq="inputs/CAMI_low/RL_S001__insert_270.fq.gz"
    output: "outputs/gs_read_mapping/{source_genome}-{contig}.fq"
    conda: "envs/seqtk.yml"
    resources:
        mem_mb = 2000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    seqtk subseq {input.fq} {input.lst} > {output}
    '''

rule split_source_genomes_by_contig:
    input: "inputs/CAMI_low/source_genomes/{source_genome}.fna"
    output: "outputs/source_genome_contigs/{source_genome}-{contig}.fa"
    resources: 
        mem_mb = 2000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    grep -A 1 {wildcards.contig} {input} > {output}
    '''

rule index_source_genome_contigs:
    input: "outputs/source_genome_contigs/{source_genome}-{contig}.fa"
    output: "outputs/bowtie2_index_source_genome_contigs/{source_genome}-{contig}.1.bt2"
    resources: 
        mem_mb = 2000,
        tmpdir=TMPDIR
    params: prefix = lambda wildcards: "outputs/bowtie2_index_source_genome_contigs/" + wildcards.source_genome + "-" + wildcards.contig
    conda: "envs/bowtie2.yml"
    threads: 1
    shell:'''
    bowtie2-build {input} {params.prefix}
    '''

rule map_reads_to_source_genome_contig:
    input:
        index="outputs/bowtie2_index_source_genome_contigs/{source_genome}-{contig}.1.bt2",
        fq = "outputs/gs_read_mapping/{source_genome}-{contig}.fq",
    output: "outputs/bowtie2_source_genome_contigs/{source_genome}-{contig}.bam"
    params: prefix = lambda wildcards: "outputs/bowtie2_index_source_genome_contigs/" + wildcards.source_genome + "-" + wildcards.contig
    resources: 
        mem_mb = 2000,
        tmpdir=TMPDIR
    conda: "envs/bowtie2.yml"
    threads: 1
    shell:'''
    bowtie2 -x {params.prefix} --interleaved {input.fq} -p {threads} --end-to-end | \
    samtools view -Sbh --threads {threads} - > {output}
    '''

rule combine_annotations_with_read_mapping_info:
    input:
        bam="outputs/bowtie2_source_genome_contigs/{source_genome}-{contig}.bam",
        gff="outputs/bakta_source_genomes/{source_genome}_no_fasta.gff"
    output:"outputs/gs_read_annotations/{source_genome}-{contig}.tsv"
    resources: 
        mem_mb = 4000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    threads: 1
    shell:'''
    intersectBed -a {input.gff} -b {input.bam} -s -bed -wa -wb > {output}
    '''

rule combine_annotated_reads_per_genome:
    input: tsv=expand("outputs/gs_read_annotations/{source_genome_and_contig}.tsv", source_genome_and_contig = SOURCE_GENOME_AND_CONTIG)
    output: tsv="outputs/outputs/gs_read_annotations/CAMI_low_gs_read_annotations.tsv"
    resources: 
        mem_mb = 256000,
        tmpdir=TMPDIR
    conda: "envs/tidyverse.yml"
    threads: 1
    script: "scripts/combine_read_mapping_annotations.R"

#rule combine_annotated_reads_per_genome_with_eggnog:

