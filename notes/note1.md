# detecting recombinant phages

experiment explanation

reference phages:

EC2D2: https://www.ncbi.nlm.nih.gov/nuccore/MZ501100.1?report=fasta

EM11: https://www.ncbi.nlm.nih.gov/nuccore/MZ501111.1?report=fasta

EM60: https://www.ncbi.nlm.nih.gov/nuccore/MZ501093.1?report=fasta

The phage reference genomes were rotated to have the terminal repeat region at the end.

## clones

Single phage isolates were taken from a recombinant population and they were sequenced. We assembled the sequencing reads to obtain the genome of the phage clone.

### SNPs distribution

We want to map on the clone genome the three references and see how the SNPs are distributed.

#### Pipeline

1. Map the ancestral sequences on the clone

    we starded considering all 3 phages but, since phage EC2D2 was never mapping we removed it from the analysis

2. Extract the SNPs from the alignment

3. Plot the SNP density of each ancestral sequence on the clone

4. summarise mismatches on a line:

    a. convolution

    b. normalization

    c. keeping values 1

    d. result