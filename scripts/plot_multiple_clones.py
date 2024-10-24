import numpy as np
import csv
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
from array_compression import decompress_array, retrive_compressed_array_from_str, compress_array
from coordinate_converter import get_hyb_ref_map


def get_mismatch_arrays(mismatch_array_path):
    csv.field_size_limit(sys.maxsize)

    ancestral_names = []
    mismatch_arrays = []
    mapping_starts = []
    mapping_ends = []

    with open(mismatch_array_path) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            # extracting mismatch arrays and mapping info
            ancestral_name = line[0]
            mapping_start = int(line[1])
            mapping_end = int(line[2])

            compressed_mismatch_array = retrive_compressed_array_from_str(line[3])
            mismatch_array = decompress_array(compressed_mismatch_array)

            ancestral_names.append(ancestral_name)
            mismatch_arrays.append(mismatch_array)
            mapping_starts.append(mapping_start)
            mapping_ends.append(mapping_end)

    return ancestral_names, mismatch_arrays, mapping_starts, mapping_ends


def npz_extract(npz_file):
    npz = np.load(npz_file)
    lst = npz.files
    for item in lst:
        array = npz[item]
    return array


def get_len_seq(fasta_path):
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return int(len(record.seq))


def fill(array, start, end, mismatch_array):
    """
    fill the array with the mismatch information starting from index "start"
    """
    for i in range(len(array)):
        if i >= start and i <= end:
            array[i] += mismatch_array[i - start]
    # return array[:start]+mismatch_array+array[start+len(mismatch_array):]
    return array


def get_mismatch_lines(
    len_seq,
    ancestral_names,
    mismatch_arrays,
    mapping_starts,
    mapping_ends,
    k,
):
    """
    create an array that summarises the normalised mismatch distribution of all the ancestral sequences
    """

    # create a dictionary with the ancestral names as keys (just one for each ancestral sequence) and the empy mismatch arrays as values,
    lines = defaultdict(None)
    for i in range(len(ancestral_names)):
        lines[ancestral_names[i]] = np.zeros(len_seq)

    # put all mismatches in one array, useful for normalization
    total_mismatches = np.zeros(len_seq)
    for i in range(len(ancestral_names)):
        total_mismatches = fill(total_mismatches, mapping_starts[i], mapping_ends[i], mismatch_arrays[i])

    # extract the mismatches of each ancestral sequence and normalise them
    for name, array in lines.items():
        for i in range(len(ancestral_names)):
            if ancestral_names[i] == name:
                array = fill(array, mapping_starts[i], mapping_ends[i], mismatch_arrays[i])

        convolved_mismatch_array = np.convolve(array, np.ones(k), "valid") / k
        convolved_total_array = np.convolve(total_mismatches, np.ones(k), "valid") / k
        normalised_array = np.divide(
            convolved_mismatch_array,
            convolved_total_array,
            out=np.zeros_like(convolved_mismatch_array),
            where=convolved_total_array != 0,
        )
        array = normalised_array  # lines now has normalised arrays for each ancestral sequence

    return lines


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--populations", help="population names")
    parser.add_argument("--clones", help="clone names")
    parser.add_argument("--hybrid_ref", help="hybrid reference")
    parser.add_argument("--k", help="convolution window")
    parser.add_argument("--out", help="output path of the plot")

    args = parser.parse_args()
    populations = args.populations.split(",")[:-1]
    clones = args.clones.split(",")[:-1]
    hybrid_ref_path = args.hybrid_ref
    k = int(args.k)
    output_path = args.out

    n_plots = len(populations) * len(clones) + 3
    plot_index = 0
    fig, axs = plt.subplots(n_plots, 1, figsize=(18, 6))

    for population in populations:
        for clone in clones:

            mismatch_array_path = f"results/mismatch_arrays/{population}/{population}_{clone}.tsv"
            # recombination_distribution_path=f'results/genomewide_recombination/{population}/{clone}.npz'

            ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)

            # recombination_distribution = npz_extract(recombination_distribution_path)

            clone_to_hybrid_alignment = f"results/alignments/evidences/{population}/{clone}/{population}_{clone}.bam"

            clone_genome_path = f"data/clones_genomes/{population}/{population}_{clone}.fasta"

            clone_len = get_len_seq(clone_genome_path)

            hyb_len = get_len_seq(hybrid_ref_path)

            mismatch_lines = get_mismatch_lines(
                clone_len, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, k
            )

            # create a single array that simplifies the normalised mismatch information (0 if the normalization is less then one, 1 or 2 in the other case depending on the ancestral reference)
            single_array_normalised_mismatches = np.zeros(clone_len)
            c = 0
            for name, array in mismatch_lines.items():
                # summarise the mismatch density of all phages
                c += 1  # zero for no evidence, 1 for EM11, 2 for EM60 #this order is determined by the order of the lines in the dictionary
                for i in range(len(array)):
                    if array[i] == 1:
                        single_array_normalised_mismatches[i] = c

            # convert to the hybrid reference coordinates
            conv_map = get_hyb_ref_map(clone_to_hybrid_alignment)
            converted_mismatch_line = np.zeros(hyb_len)
            for clone_coord in range(len(single_array_normalised_mismatches)):
                hyb_coord = conv_map[clone_coord]
                if hyb_coord == None:
                    continue  # i have to create a mask
                converted_mismatch_line[hyb_coord] = single_array_normalised_mismatches[clone_coord]

            # plot function
            axs[plot_index].plot(converted_mismatch_line, label=f"{population} {clone}")
            axs[plot_index].set_title(f"{population} {clone}")
            axs[plot_index].legend()
            plot_index += 1

    plt.savefig(output_path, bbox_inches="tight")
