configfile: 'config.yaml'

pars=config['Parameters']
depths_list=config['Depths']
individuals_list=config['Individuals']
input_paths=config['Input_Paths']
n_individuals = len(individuals_list)

def get_triangular_bcf_comparisons():
    triangular_pair_list = []
    for a in range(1, n_individuals):
        for b in range(a+1, n_individuals+1):
            triangular_pair_list.append({
                "indi_a": a,
                "indi_b": b
            })
    return triangular_pair_list

def get_out():
    trees=expand("RAxML_trees/aligned_{depth}X.raxml.bestTree", depth = depths_list)
    return [trees,"distAngsd_results/distAngsd_results.txt"]

rule all:
    input: get_out()

#simulate sequencing with art
rule art_generate_reads:
    input:
        fasta = "true_fastas/{individual}_{haplotype}.fasta"
    params:
        l=pars['art_l'],
        m=pars['art_m'],
        s=pars['art_s'],
        prefix="haploid_fastas/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X"
    output:
        haploid_R1="haploid_fastas/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X_R1.fq",
        haploid_R2="haploid_fastas/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X_R2.fq"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate diptest

        half_depth=$(echo "{wildcards.depth} / 2" | bc -l)
        art_illumina \
            -ss HS25 \
            -i {input.fasta} \
            -l {params.l} \
            -f "$half_depth" \
            -p \
            -m {params.m} \
            -s {params.s} \
            -rs 374 \
            -o {params.prefix}

        mv {params.prefix}1.fq {params.prefix}_R1.fq
        mv {params.prefix}2.fq {params.prefix}_R2.fq
        echo $half_depth
        '''

#pair haplotypes per individual
rule pair_diploids:
    input: 
        hap1_read1 = "haploid_fastas/haploid_depth_{depth}X/individual_{individual}_hap1_depth_{depth}X_R1.fq",
        hap1_read2 = "haploid_fastas/haploid_depth_{depth}X/individual_{individual}_hap1_depth_{depth}X_R2.fq",
        hap2_read1 = "haploid_fastas/haploid_depth_{depth}X/individual_{individual}_hap2_depth_{depth}X_R1.fq",
        hap2_read2 = "haploid_fastas/haploid_depth_{depth}X/individual_{individual}_hap2_depth_{depth}X_R2.fq",
    output: 
        dipout_read1 = "diploid_fastas/depth_{depth}X/individual_{individual}_depth_{depth}X_R1.fq",
        dipout_read2 = "diploid_fastas/depth_{depth}X/individual_{individual}_depth_{depth}X_R2.fq"
    shell:
        '''
        cat {input.hap1_read1} {input.hap2_read1} > {output.dipout_read1}
        cat {input.hap1_read2} {input.hap2_read2} > {output.dipout_read2}
        '''

#alignment with bwa mem
rule bwa_generate_bams:
    input: 
        read1 = "diploid_fastas/depth_{depth}X/individual_{individual}_depth_{depth}X_R1.fq",
        read2 = "diploid_fastas/depth_{depth}X/individual_{individual}_depth_{depth}X_R2.fq"
    params:
        ref=input_paths['reference']
    threads:2
    output:
        bam = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate alignenv

        bwa mem -t{threads}  {params.ref} {input.read1} {input.read2} |
        samtools view -b - |
        samtools sort -o {output.bam}
        samtools index {output.bam}
        '''

#mpileup with bcftools
rule bcftools_mpileup:
    input: 
        bam_a="bams/depth_{depth}X/individual_{indi_a}_depth_{depth}X.sorted.bam",
        bam_b="bams/depth_{depth}X/individual_{indi_b}_depth_{depth}X.sorted.bam"
        
    params:
        min_MQ=pars['bcftools_min_MQ'],
        min_BQ=pars['bcftools_min_BQ'],
        ref=input_paths['reference']
    threads: 3
    output: 
        output_bcf="mpileups/depth_{depth}X/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}X.bcf",
        output_vcf_gz="mpileups/depth_{depth}X/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}X.vcf.gz"
    shell:
        '''
        module load perl/5.38.0
        module load gsl/2.5 
        module load bcftools/1.21

        bcftools mpileup -Ou --min-MQ {params.min_MQ} --min-BQ {params.min_BQ} -f {params.ref} {input.bam_a} {input.bam_b} --threads {threads} | 
        bcftools filter -Ou -e 'INFO/INDEL!=0' -o {output.output_bcf}
        bcftools view {output.output_bcf} -Oz -o {output.output_vcf_gz}
        '''

#calculate distance with distAngsd
rule distAngsd_distance_calc:
    input:
        vcf="mpileups/depth_{depth}X/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}X.vcf.gz"
    params:
        distAngsd_path = input_paths['distAngsd_path']
    output:
        tmp_out="distAngsd_results/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}.txt"
    shell:
        '''
        module load lapack
        estimate=$({params.distAngsd_path}/distAngsd -vcf {input.vcf} -model GTR -par 1.0,2.0,1.0,1.0,2.0,1.0,0.25,0.25,0.25 2>&1 | grep "Estimated t =" | awk '{{print $4}}')
        echo -e "{wildcards.depth}\t{wildcards.indi_a}\t{wildcards.indi_b}\t${{estimate}}" > {output.tmp_out}
        '''

rule build_distAngsd_results_file:
    input: 
        all_tmp_outs = expand([f"distAngsd_results/individuala_{pair['indi_a']}_individualb_{pair['indi_b']}_depth_{{depth}}.txt" for pair in get_triangular_bcf_comparisons()], depth=depths_list)
    output:
        file="distAngsd_results/distAngsd_results.txt"
    run:
        with open(output.file, 'w') as out_f:
            out_f.write("depth\tindi_a\tindi_b\tdistance\n")
            for tmp_out in input.all_tmp_outs:
                with open(tmp_out) as f:
                    out_f.write(f.read())

#create fastas with angsd
rule angsd_generate_fastas:
    input:
        bam = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam"
    params:
        min_MQ=pars['angsd_min_MQ'],
        min_BQ=pars['angsd_min_BQ'],
        ref=input_paths['reference'],
        proper_pairs=pars['angsd_only_proper_pairs'],
        prefix="angsd_fastas/depth_{depth}/individual_{individual}_depth_{depth}"

    threads: 3
    output:
        fasta="angsd_fastas/depth_{depth}/individual_{individual}_depth_{depth}.fa"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate angsd

        angsd -i {input.bam} -ref {params.ref} -out {params.prefix} -nThreads {threads} -doFasta 2 -doCounts 1 -minQ {params.min_BQ} -minMapQ {params.min_MQ} -only_proper_pairs {params.proper_pairs}
        gunzip "{params.prefix}.fa.gz"
        filename=$(basename {output.fasta} .fa)
        id=$(echo "$filename" | cut -d'_' -f2)
        tmp_file="{output.fasta}.tmp"
        awk -v p=">$id" '/^>/ {{ print p; next }} {{ print }}' {output.fasta} > "$tmp_file"
        mv "$tmp_file" {output.fasta}
        '''
#concat fastas
rule concat_fastas:
    input:
        fasta_list=expand("angsd_fastas/depth_{{depth}}/individual_{individual}_depth_{{depth}}.fa", individual = individuals_list)
    params:
        ref=input_paths['reference'],
    output:
        aligned_fasta="aligned_fastas/aligned_depth_{depth}X.fasta"
    run:
        with open(output.aligned_fasta, 'w') as out_f:
            for indi_fasta in input.fasta_list:
                with open(indi_fasta) as f:
                    out_f.write(f.read())
            with open(params.ref) as ref_f:
                out_f.write(ref_f.read())

#generate trees with RAxML
rule RAxML_tree_gen:
    input:
        aligned_fasta = "aligned_fastas/aligned_depth_{depth}X.fasta"
    params:
        prefix="RAxML_trees/aligned_{depth}X",
        outgroup_name=input_paths['outgroup_name']
    threads: pars['RAxML_threads']
    output:
        "RAxML_trees/aligned_{depth}X.raxml.bestTree"

    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate RAxML

        raxml-ng \
        --all \
        --tree rand{{50}},pars{{50}} \
        --seed 12345 \
        --msa {input.aligned_fasta} \
        --model 'GTR{{1.0,2.0,1.0,1.0,2.0,1.0}}+FE' \
        --threads {threads} \
        --prefix {params.prefix} \
        --outgroup {params.outgroup_name}
        '''