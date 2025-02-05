# recombinant_phages_study

This repository contains a pipeline used to analyse recombinant phage clones.

The plots generated are described [here](documentation/plots.md).

All the steps of the pipeline and the functioning of the scripts are described [here](documentation/scripts.md).

# Configuration

The [config.yml](/config.yml) file contains all the parameters of the pipeline:

- HMM: set the cores available for the HMM algorithm and the parameters of the HMM model used to create the plots

- Optimization of recombination parameter: the optimization of the transition parameter is run independently from the whole pipeline, it consists of several runs of the HMM algorithm over all the clones. At each iteration a different value of transition probability is used and the total log likelihood is saved. The output of the optimization is just a plot that shows the average trend of log likelihood depending on the recombination parameter. The flag parameter is a boolean value that defines if the optimization will be executed or not. The parameter "values" takes comma separated entries for the transition probability.

- Plots: these parameters control the final plot of the pipeline, containing all the clones aligned to the hybrid reference.
    - convolution_window: convolution window applied to the mismatch distribution, affects the resolution of the clone genome line.
    - common_x_axis: not important, probably will be removed
    - interline: space in between lines representing clone genomes
    - thickness: thickness of the line representing the clone genome
    - colors:

## Input folder

The input folder should have the following structure:

<pre>
data/
    clones_genomes/
        [population_code]/
            [population_code]_[clone_id].fasta
            ...
    references.fasta
</pre>

# Running the pipeline

## Local execution

<pre>
snakemake --profile local all
</pre>

## HPC execution 

<pre>
snakemake --profile cluster all
</pre>
