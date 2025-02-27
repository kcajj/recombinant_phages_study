population_names=["P1","P2","P3"]
clone_names=["C1","C2","C3","C4"]

populations_string = ""
for population in population_names:
    populations_string += population + ","
clones_string = ""
for clone in clone_names:
    clones_string += clone + ","

configfile: "config.yml"

hmm = config["hmm"]
x_axis = config["plots"]["common_x_axis"]
k = config["plots"]["convolution_window"]
mismatch_threshold = config["plots"]["mismatch_threshold"]
interline = config["plots"]["interline"]
thickness = config["plots"]["thickness"]
colors = config["plots"]["colors"]

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
        plots='results/plots/mismatch_density/{population}_{isolate}.pdf',
    conda:
        'conda_envs/sci_py.yml'
    params:
        x_axis = x_axis,
        k = k
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
        plots='results/plots/mismatch_line/{population}_{isolate}.pdf',
    conda:
        'conda_envs/sci_py.yml'
    params:
        x_axis = x_axis,
        k = k,
        mismatch_threshold = mismatch_threshold,
    shell:
        """
        python scripts/plot_mismatch_line.py \
            --clone {input.clone} \
            --mismatches {input.mismatch} \
            --x_axis {params.x_axis} \
            --k {params.k} \
            --mismatch_threshold {params.mismatch_threshold} \
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
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/extract_evidence_arrays.py \
            --bam {input.bam} \
            --msa_refs {input.msa} \
            --evidences_out {output.evidences} \
        """

rule prediction_arrays:
    input:
        evidences = rules.evidence_arrays.output.evidences,
        hybrid_ref = rules.hybrid_ref.output.hybrid_ref
    output:
        predictions='results/prediction_arrays/{population}/{population}_{isolate}.tsv',
        coverage='results/coverage_arrays/{population}/{population}_{isolate}.npz'
    conda:
        'conda_envs/sci_py.yml'
    params:
        cores = hmm["cores"],
        initial_probability = hmm["initial_probability"]["A"]+","+hmm["initial_probability"]["B"],
        transition_probability = hmm["transition_probability"]["A"]["A"]+","+hmm["transition_probability"]["A"]["B"]+"/"+hmm["transition_probability"]["B"]["A"]+","+hmm["transition_probability"]["B"]["B"],
        emission_probability = hmm["emission_probability"]["A"][0]+","+hmm["emission_probability"]["A"][1]+","+hmm["emission_probability"]["A"][2]+"/"+hmm["emission_probability"]["B"][0]+","+hmm["emission_probability"]["B"][1]+","+hmm["emission_probability"]["B"][2]
    shell:
        """
        python scripts/hmm_prediction_arrays.py \
            --evidences {input.evidences} \
            --hybrid_ref {input.hybrid_ref} \
            --predictions_out {output.predictions} \
            --coverage_out {output.coverage} \
            --cores {params.cores} \
            --initial_p {params.initial_probability}\
            --transition_p {params.transition_probability}\
            --emission_p {params.emission_probability}
        """

rule plot_references_coverage:
    input:
        coverage = rules.prediction_arrays.output.coverage,
        msa = rules.msa.output.msa
    output:
        plots='results/plots/references_coverage/{population}/{population}_{isolate}.pdf',
    conda:
        'conda_envs/sci_py.yml'
    shell:
        """
        python scripts/plot_references_coverage.py \
            --msa_refs {input.msa} \
            --coverage {input.coverage} \
            --out {output.plots}
        """

rule clones_processing:
    input:
        mismatch_density_plots=expand(rules.plot_mismatch_density.output.plots, population=population_names, isolate=clone_names),
        mismatch_line_plots=expand(rules.plot_mismatch_line.output.plots, population=population_names, isolate=clone_names),
        coverage_plots=expand(rules.plot_references_coverage.output.plots, population=population_names, isolate=clone_names),
    output:
        finish = 'results/all_clones_processed.txt'
    shell:
        """
        touch {output.finish}
        """

rule map_ancestral_to_hybridref:
    input:
        ancestral = ancestral_phages,
        hybrid = rules.hybrid_ref.output.hybrid_ref
    output:
        sam = 'results/alignments/ancestral_to_hybridref/ancestral_to_hybridref.sam',
        bam = 'results/alignments/ancestral_to_hybridref/ancestral_to_hybridref.bam',
        bai = 'results/alignments/ancestral_to_hybridref/ancestral_to_hybridref.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a {input.hybrid} {input.ancestral} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """

rule plot_multiple_clones:
    input:
        finish = rules.clones_processing.output.finish,
        hybrid_ref = rules.hybrid_ref.output.hybrid_ref,
        ancestral_alignment = rules.map_ancestral_to_hybridref.output.bam,
    output:
        plot='results/plots/multiple_clones_plot.pdf',
    conda:
        'conda_envs/sci_py.yml'
    params:
        populations = populations_string,
        clones = clones_string,
        k = k,
        mismatch_threshold = mismatch_threshold,
        interline = interline,
        thickness = thickness,
        colors = colors
    shell:
        """
        python scripts/plot_multiple_clones.py \
            --populations {params.populations} \
            --clones {params.clones} \
            --hybrid_ref {input.hybrid_ref} \
            --ancestral_alignment {input.ancestral_alignment} \
            --k {params.k} \
            --mismatch_threshold {params.mismatch_threshold} \
            --interline {params.interline} \
            --thickness {params.thickness} \
            --colors {params.colors} \
            --out {output.plot}
        """

rule optimize_recombination_parameter:
    input:
        finish = rules.clones_processing.output.finish,
    output:
        plot = 'results/plots/parameter_optimization.pdf'
    conda:
        'conda_envs/sci_py.yml'
    params:
        populations = populations_string,
        clones = clones_string,
        cores = hmm["cores"],
        initial_probability = hmm["initial_probability"]["A"]+","+hmm["initial_probability"]["B"],
        transition_probability = config["optimization_recombination_parameter"]["values"],
        emission_probability = hmm["emission_probability"]["A"][0]+","+hmm["emission_probability"]["A"][1]+","+hmm["emission_probability"]["A"][2]+"/"+hmm["emission_probability"]["B"][0]+","+hmm["emission_probability"]["B"][1]+","+hmm["emission_probability"]["B"][2]
    shell:
        """
        python scripts/optimize_recombination_parameter.py \
            --populations {params.populations} \
            --clones {params.clones} \
            --out {output.plot} \
            --cores {params.cores} \
            --initial_p {params.initial_probability}\
            --transition_p {params.transition_probability}\
            --emission_p {params.emission_probability}
        """

optimize = []
if config["optimization_recombination_parameter"]["flag"]:
    optimize.append(rules.optimize_recombination_parameter.output.plot)

rule all:
    input:
        finish = rules.clones_processing.output.finish,
        plot = rules.plot_multiple_clones.output.plot,
        optimization = optimize
    shell:
        """
        rm {input.finish}
        """