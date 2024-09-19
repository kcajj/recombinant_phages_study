clones = 'data/clones_genomes/{population}/{population}_{isolate}.fasta'
ancestral_phages = 'data/references.fasta'

rule map_ancestral_to_clones:
    input:
        clone = clones,
        ancestral = ancestral_phages
    output:
        sam = 'results/clones/alignments/{population}/{population}_{isolate}.sam',
        bam = 'results/clones/alignments/{population}/{population}_{isolate}.bam',
        bai = 'results/clones/alignments/{population}/{population}_{isolate}.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a {input.clone} {input.ancestral} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """

rule mismatch_arrays:
    input:
        bam = rules.map_ancestral_to_clones.output.bam,
        clone = clones
    output:
        mismatch='results/clones/mismatch_arrays/{population}/{population}_{isolate}.tsv',
        coverage='results/clones/coverage_arrays/{population}/{population}_{isolate}.npz'
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/extract_mismatch_arrays.py \
            --bam {input.bam} \
            --clone {input.clone} \
            --mismatch_out {output.mismatch} \
            --coverage_out {output.coverage}
        """

rule all:
    input:
        mismatch=expand(rules.mismatch_arrays.output.mismatch, population=["P1"], isolate=["C1"])