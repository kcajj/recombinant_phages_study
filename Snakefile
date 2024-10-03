configfile: "config.yml"

clones = 'data/clones_genomes/{population}/{population}_{isolate}.fasta'
ancestral_phages = 'data/references.fasta'

rule map_ancestral_to_clones:
    input:
        clone = clones,
        ancestral = ancestral_phages
    output:
        sam = 'results/alignments/mismatches/{population}/{population}_{isolate}.sam',
        bam = 'results/alignments/mismatches/{population}/{population}_{isolate}.bam',
        bai = 'results/alignments/mismatches/{population}/{population}_{isolate}.bam.bai'
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
        mismatch='results/mismatch_arrays/{population}/{population}_{isolate}.tsv'
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/extract_mismatch_arrays.py \
            --bam {input.bam} \
            --clone {input.clone} \
            --mismatch_out {output.mismatch}
        """

rule plot_mismatch_density:
    input:
        clone = clones,
        mismatch = rules.mismatch_arrays.output.mismatch
    output:
        plots='results/plots/mismatch_density/{population}_{isolate}.png',
    conda:
        'conda_envs/sci_py.yml'
    params:
        x_axis = config["common_x_axis"],
        k = config["convolution_window"]
    shell:
        """
        python scripts/plot_mismatch_density.py \
            --clone {input.clone} \
            --mismatches {input.mismatch} \
            --x_axis {params.x_axis} \
            --k {params.k} \
            --out {output.plots}
        """

rule plot_mismatch_line:
    input:
        clone = clones,
        mismatch = rules.mismatch_arrays.output.mismatch
    output:
        plots='results/plots/mismatch_line/{population}_{isolate}.png',
    conda:
        'conda_envs/sci_py.yml'
    params:
        x_axis = config["common_x_axis"],
        k = config["convolution_window"]
    shell:
        """
        python scripts/plot_mismatch_line.py \
            --clone {input.clone} \
            --mismatches {input.mismatch} \
            --x_axis {params.x_axis} \
            --k {params.k} \
            --out {output.plots}
        """

rule msa:
    input:
        reference = ancestral_phages
    output:
        msa = 'results/msa/msa_refs.fasta'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        mafft --auto \
            {input.reference} \
            > {output.msa}
        """

rule hybrid_ref:
    input:
        msa = rules.msa.output.msa
    output:
        hybrid_ref = 'results/msa/hybrid_ref.fasta'
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/hybrid_reference.py \
            --msa {input.msa} \
            --out {output.hybrid_ref}
        """

rule mapping:
    input:
        clone=clones,
        ref=rules.hybrid_ref.output.hybrid_ref
    output:
        sam = 'results/alignments/evidences/{population}/{isolate}/{population}_{isolate}.sam',
        bam = 'results/alignments/evidences/{population}/{isolate}/{population}_{isolate}.bam',
        bai = 'results/alignments/evidences/{population}/{isolate}/{population}_{isolate}.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a {input.ref} {input.clone} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """

rule evidence_arrays:
    input:
        bam = rules.mapping.output.bam,
        msa = rules.msa.output.msa
    output:
        evidences='results/evidence_arrays/{population}/{population}_{isolate}.tsv',
        coverage='results/coverage_arrays/{population}/{population}_{isolate}.npz'
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/extract_evidence_arrays.py \
            --bam {input.bam} \
            --msa_refs {input.msa} \
            --evidences_out {output.evidences} \
            --coverage_out {output.coverage}
        """

rule prediction_arrays:
    input:
        evidences = rules.evidence_arrays.output.evidences,
        msa = rules.msa.output.msa
    output:
        predictions='results/prediction_arrays/{population}/{population}_{isolate}.tsv'
    conda:
        'conda_envs/sci_py.yml'
    params:
        cores = config["cores"],
        initial_probability = config["initial_probability"]["A"]+","+config["initial_probability"]["B"],
        transition_probability = config["transition_probability"]["A"]["A"]+","+config["transition_probability"]["A"]["B"]+"/"+config["transition_probability"]["B"]["A"]+","+config["transition_probability"]["B"]["B"],
        emission_probability = config["emission_probability"]["A"][0]+","+config["emission_probability"]["A"][1]+","+config["emission_probability"]["A"][2]+"/"+config["emission_probability"]["B"][0]+","+config["emission_probability"]["B"][1]+","+config["emission_probability"]["B"][2]
    shell:
        """
        python scripts/hmm_prediction_arrays.py \
            --evidences {input.evidences} \
            --msa_refs {input.msa} \
            --out {output.predictions} \
            --cores {params.cores} \
            --initial_p {params.initial_probability}\
            --transition_p {params.transition_probability}\
            --emission_p {params.emission_probability}
        """

rule plot_references_coverage:
    input:
        predictions = rules.prediction_arrays.output.predictions,
        msa = rules.msa.output.msa
    output:
        plots='results/plots/references_coverage/{population}/{population}_{isolate}.png',
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/plot_references_coverage.py \
            --predictions {input.predictions} \
            --msa_refs {input.msa} \
            --out {output.plots}
        """

rule all:
    input:
        mismatch_density_plots=expand(rules.plot_mismatch_density.output.plots, population=["P1","P2","P3"], isolate=["C1","C2","C3","C4"]),
        mismatch_line_plots=expand(rules.plot_mismatch_line.output.plots, population=["P1","P2","P3"], isolate=["C1","C2","C3","C4"]),
        coverage_plots=expand(rules.plot_references_coverage.output.plots, population=["P1","P2","P3"], isolate=["C1","C2","C3","C4"]),