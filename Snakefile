configfile: 'config.yaml'

from datetime import datetime
current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

pars=config['Parameters']
depths_list=config['Depths']
input_paths=config['Input_Paths']
n_individuals=pars['n_individuals']
run_id=pars['run_id']
individuals_list=list(range(1,n_individuals+1))

def get_triangular_bcf_comparisons():
    n_individuals=pars['n_individuals']
    triangular_pair_list = []
    for a in range(1, n_individuals):
        for b in range(a+1, n_individuals+1):
            triangular_pair_list.append({
                "indi_a": a,
                "indi_b": b
            })
    return triangular_pair_list

def get_out():
    return [
        "multiqc/multiqc_report.html", 
        "flagstats_summary/flagstat_summary.tsv" 
        #"distance_error_by_2methods.png", 
        #"individual_tree.png", 
        #"haplotype_tree.png"
        ]

rule all:
    input: get_out()

rule generate_trees:
    params:
        n_individuals=pars['n_individuals'],
        branch_length_m=pars['tree_branch_length_m'],
        branch_length_sd=pars['tree_branch_length_sd'],
        hap_distance=pars['tree_hap_distance']
    output:
        haplotype_tree="trees/haplotype_tree.nwk",
        individual_tree="trees/individual_tree.nwk"
    shell:
        '''
        module load gcc/13.2.0
        module load openjdk/20.0.0
        module load R
        Rscript generate_trees.R \
        {params.n_individuals} \
        {params.branch_length_m} \
        {params.branch_length_sd} \
        {params.hap_distance} \
        {output.haplotype_tree} \
        {output.individual_tree}
        '''

rule generate_true_fastas:
    input:
        haplotype_tree="trees/haplotype_tree.nwk"
    params:
        model=pars['true_fastas_model'],
        l=pars['true_fastas_length'],
        f=pars['true_fastas_initial_freq'],
        r=pars['true_fastas_rates'],
        phy_path="true_fastas/simulated_sequences.phy",
        fasta_path="true_fastas/simulated_sequences.fasta",
        n_individuals=n_individuals
    output:
        fasta_list=expand("true_fastas/{{individual}}_{haplotype}.fasta",haplotype=["hap1", "hap2"]),
        aligned_true_fasta = "true_fastas/simulated_sequences.fasta" 
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate diptest 

        seq-gen \
        -z12345 \
        -m{params.model} \
        -l{params.l} \
        -f{params.f} \
        -r{params.r} \
        <"{input.haplotype_tree}"> \
        "{params.phy_path}"

        seqmagick convert "{params.phy_path}" "{params.fasta_path}"
        seqkit split -i -O true_fastas "{params.fasta_path}" --by-id-prefix ""   
        '''
        
#isolate and index outgroup:
rule outgroup:
    input:
        fasta="true_fastas/out_hap1.fasta"
    output:
        outgroup="outgroup/out_hap1.fasta",
        indexed_outgroup="outgroup/out_hap1.fasta.fai"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate alignenv
        module load samtools 

        cp {input.fasta} outgroup/
        bwa index {output.outgroup}
        samtools faidx {output.outgroup}
        '''

rule heatmap:
    input: 
        aligned_true_fasta = "true_fastas/simulated_sequences.fasta"
    output:
        distance.tsv = true_fastas/fasta_snp_distance.tsv
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate snp-dists

        snp-dists -m true_fastas/simulated_sequences.fasta > true_fastas/fasta_snp_distance.tsv
        '''

        
#simulate sequencing with ngsngs
rule ngsngs_generate_reads:
    input:
        fasta = "true_fastas/{individual}_{haplotype}.fasta"
    params:
        ngsngs_path = input_paths['ngsngs_path'],
        q1_path = input_paths['ngsngs_q1_path'],
        q2_path = input_paths['ngsngs_q2_path'],
        l = pars['ngsngs_fragment_legnth'],
        cl = pars['ngsngs_read_length'],
        prefix="haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X"
    output:
        haploid_R1="haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X_R1.fq",
        haploid_R2="haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_{haplotype}_depth_{depth}X_R2.fq"
    shell:
        '''
        half_depth=$(echo "{wildcards.depth} / 2" | bc -l)

        {params.ngsngs_path}/ngsngs \
        -l {params.l} \
        -cl {params.cl} \
        -seq "PE" \
        -f fq \
        -q1 {params.q1_path} \
        -q2 {params.q2_path} \
        -c $half_depth \
        -i {input.fasta} \
        -o {params.prefix}
        '''

#pair haplotypes per individual
rule pair_diploids:
    input: 
        hap1_read1 = "haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_hap1_depth_{depth}X_R1.fq",
        hap1_read2 = "haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_hap1_depth_{depth}X_R2.fq",
        hap2_read1 = "haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_hap2_depth_{depth}X_R1.fq",
        hap2_read2 = "haploid_fastqs/haploid_depth_{depth}X/individual_{individual}_hap2_depth_{depth}X_R2.fq"
    params:
        collapsed_read_prefix = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X"
    output: 
        dipout_read1 = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X_R1.fq",
        dipout_read2 = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X_R2.fq",
        collapsed_read = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X.collapsed"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate AdapterRemoval

        cat {input.hap1_read1} {input.hap2_read1} > {output.dipout_read1}
        cat {input.hap1_read2} {input.hap2_read2} > {output.dipout_read2}

        AdapterRemoval \
        --file1 {output.dipout_read1} \
        --file2 {output.dipout_read2} \
        --basename {params.collapsed_read_prefix} \
        --collapse
        '''

#check sequencing
rule fastqc:
    input:
        collapsed_read = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X.collapsed"
    output:
        fastqcs = "fastqcs/individual_{individual}_depth_{depth}X.collapsed_fastqc.html"
    shell:
        '''
        module load --auto fastqc
        fastqc {input.collapsed_read} --outdir fastqcs
        '''
rule multiqc:
    input:
        fastqcs = expand("fastqcs/individual_{individual}_depth_{depth}X.collapsed_fastqc.html", depth=depths_list,individual = individuals_list)
    output:
        multiqc = "multiqc/multiqc_report.html"
    shell:
        '''
        module load --auto multiqc
        multiqc fastqcs --outdir multiqc -f
        '''

#alignment with bwa aln
rule bwa_generate_bams:
    input: 
        collapsed_read = "diploid_fastqs/depth_{depth}X/individual_{individual}_depth_{depth}X.collapsed",
        ref="outgroup/out_hap1.fasta",
        indexed_outgroup="outgroup/out_hap1.fasta.fai"
    threads:2
    output:
        bam = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam",
        bai = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam.bai"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate alignenv

        bwa mem -t{threads} {input.ref} {input.collapsed_read} | \
        samtools view -b | samtools sort -o {output.bam} 

        samtools index {output.bam}
        '''

#check mapping
rule individual_flagstats:
    input: 
        bam = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam"
    output:
        tmp_out = temp("flagstats/individual_{individual}_depth_{depth}X_tmp.txt"),
        individual_flagstats = "flagstats/individual_{individual}_depth_{depth}X_flagstat.txt"
    shell:
        '''
        module load --auto samtools
        
        samtools flagstat {input.bam} > {output.tmp_out}
        total=$(grep "in total" "{output.tmp_out}" | awk '{{print $1}}')
        mapped=$(grep "mapped (" -m 1 "{output.tmp_out}" | awk '{{print $1}}')
        mapped_pct=$(grep "mapped (" -m 1 "{output.tmp_out}" | awk -F'[()%]' '{{print $2}}')
        echo -e "{wildcards.individual}\t{wildcards.depth}\t$total\t$mapped\t$mapped_pct" > {output.individual_flagstats}
        '''

rule compile_flagstats:
    input:
        flagstats_list = expand("flagstats/individual_{individual}_depth_{depth}X_flagstat.txt", individual=individuals_list,depth=depths_list)
    output:
        file = "flagstats_summary/flagstat_summary.tsv"
    run:
        with open(output.file, 'w') as out_f:
            out_f.write("individual\tdepth\ttotal\tmapped\tmapped_pct\n")
            for tmp_out in input.flagstats_list:
                with open(tmp_out) as f:
                    out_f.write(f.read())


#mpileup with bcftools
rule bcftools_mpileup:
    input: 
        bam_a="bams/depth_{depth}X/individual_{indi_a}_depth_{depth}X.sorted.bam",
        bam_b="bams/depth_{depth}X/individual_{indi_b}_depth_{depth}X.sorted.bam"
        
    params:
        min_MQ=pars['bcftools_min_MQ'],
        min_BQ=pars['bcftools_min_BQ'],
        ref="outgroup/out_hap1.fasta"
    output: 
        output_bcf="mpileups/depth_{depth}X/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}X.bcf",
        output_vcf_gz="mpileups/depth_{depth}X/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}X.vcf.gz"
    shell:
        '''
        module load perl/5.38.0
        module load gsl/2.5 
        module load bcftools/1.21

        bcftools mpileup -Ou --min-MQ {params.min_MQ} --min-BQ {params.min_BQ} -f {params.ref} {input.bam_a} {input.bam_b} | 
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
        tmp_out="distAngsd_results/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}.txt",
        tmp_log="distAngsd_results/individuala_{indi_a}_individualb_{indi_b}_depth_{depth}_log.txt"
    shell:
        '''
        module load lapack

        {params.distAngsd_path}/distAngsd \
        -vcf {input.vcf} \
        -model GTR \
        -par 1.0,2.0,1.0,1.0,2.0,0.25,0.25,0.25,0.25 > {output.tmp_log}

        estimate=$(grep "Estimated t =" {output.tmp_log} | awk '{{print $4}}')
        echo -e "{wildcards.depth}\t{wildcards.indi_a}\t{wildcards.indi_b}\t${{estimate}}" > {output.tmp_out}
        '''

rule build_distAngsd_results_file:
    input: 
        all_tmp_outs = expand([f"distAngsd_results/individuala_{pair['indi_a']}_individualb_{pair['indi_b']}_depth_{{depth}}.txt" for pair in get_triangular_bcf_comparisons()], depth=depths_list),
        all_logs = expand([f"distAngsd_results/individuala_{pair['indi_a']}_individualb_{pair['indi_b']}_depth_{{depth}}_log.txt" for pair in get_triangular_bcf_comparisons()], depth=depths_list)
    output:
        file="distAngsd_results/distAngsd_results.txt",
        log_file="distAngsd_results/all_logs.txt"
    run:
        with open(output.file, 'w') as out_f:
            out_f.write("depth\tindi_a\tindi_b\tdistance\n")
            for tmp_out in input.all_tmp_outs:
                with open(tmp_out) as f:
                    out_f.write(f.read())
        
        with open(output.log_file, 'w') as out_f:
            out_f.write("depth\tindi_a\tindi_b\tdistance\n")
            for log_file in input.all_logs:
                with open(log_file) as f:
                    out_f.write(f.read())
        

#create fastas with angsd
rule angsd_generate_fastas:
    input:
        bam = "bams/depth_{depth}X/individual_{individual}_depth_{depth}X.sorted.bam"
    params:
        min_MQ=pars['angsd_min_MQ'],
        min_BQ=pars['angsd_min_BQ'],
        ref="outgroup/out_hap1.fasta",
        prefix="angsd_fastas/depth_{depth}/individual_{individual}_depth_{depth}"

    threads: 3
    output:
        fasta="angsd_fastas/depth_{depth}/individual_{individual}_depth_{depth}.fa"
    shell:
        '''
        module load miniconda
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate angsd

        angsd \
        -i {input.bam} \
         -ref {params.ref} \
         -out {params.prefix} \
         -nThreads {threads} \
         -doFasta 2 \
         -doCounts 1 \
         -minQ {params.min_BQ} \
         -minMapQ {params.min_MQ} 
        
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
        ref="outgroup/out_hap1.fasta",
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
        outgroup_name='out_hap1'
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

rule plot_results:
    input: 
        trees=expand("RAxML_trees/aligned_{depth}X.raxml.bestTree", depth = depths_list),
        distAngsd_resultsfile="distAngsd_results/distAngsd_results.txt",
        individual_tree="trees/individual_tree.nwk",
        haplotype_tree="trees/haplotype_tree.nwk"
    output:
        "plots/distance_error_by_2methods.png",
        "plots/individual_tree.png",
        "plots/haplotype_tree.png"
    shell:
        '''
        module load gcc/13.2.0
        module load openjdk/20.0.0
        module load R
        Rscript plotting.R {current_time} {run_id}
        '''

        