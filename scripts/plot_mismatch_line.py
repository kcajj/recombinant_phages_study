import numpy as np
import csv
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
from array_compression import decompress_array, retrive_compressed_array_from_str, compress_array


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


def summarise_normalise_mismatch_arrays(len_seq, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, k):
    """
    create an array that summarises the normalised mismatch distribution of all the ancestral sequences
    """

    # create a dictionary with the ancestral names as keys (just one for each ancestral sequence) and the empy mismatch arrays as values,
    distributions = defaultdict(None)
    for i in range(len(ancestral_names)):
        distributions[ancestral_names[i]] = np.zeros(len_seq)

    # put all mismatches in one array, useful for normalization
    total_mismatches = np.zeros(len_seq)
    for i in range(len(ancestral_names)):
        total_mismatches = fill(total_mismatches, mapping_starts[i], mapping_ends[i], mismatch_arrays[i])

    # extract the mismatches of each ancestral sequence and normalise them
    for name, array in distributions.items():
        for i in range(len(ancestral_names)):
            if ancestral_names[i] == name:
                distributions[name] = fill(array, mapping_starts[i], mapping_ends[i], mismatch_arrays[i])
        # convolution and normalisation
        convolved_mismatch_array = np.convolve(array, np.ones(k), "valid") / k
        convolved_total_array = np.convolve(total_mismatches, np.ones(k), "valid") / k
        normalised_array = np.divide(
            convolved_mismatch_array,
            convolved_total_array,
            out=np.zeros_like(convolved_mismatch_array),
            where=convolved_total_array != 0,
        )

        distributions[name] = normalised_array  # the distributions dictionary stores the normalised arrays

    return distributions


def plot_mismatch_line(
    distributions,
    mismatch_line,
    population,
    isolate,
    k,
    x_axis,
    out_folder,
):

    # create a colour array with 0=gray and from 1 onwards C0, C1, C2, C3 ecc
    colors = {"C0": "gray"}
    c = 0
    for distribution in distributions.keys():
        # for each color CX with X>0 we assign the color names starting from C0
        colors["C" + str(c + 1)] = "C" + str(c)
        c += 1
        # cycling through the lines in the same order as we did in the previous loop, to guarantee correct colours (if a line is first it gets c0)

    # create the figure
    fig, ax = plt.subplots(nrows=1, figsize=(18, 6))

    # compress the array so that it is easy to create rectangles
    plotting_array = compress_array(mismatch_line)

    # plot the rectangles
    height = 1
    y = -0.5
    c = 0
    for lenn, typee in plotting_array:
        x = c
        width = lenn
        c += lenn

        # invert colors
        if typee == 1:
            typee = 2
        elif typee == 2:
            typee = 1

        rectangle = mpatches.Rectangle((x, y), width, height, color=colors["C" + str(typee)])
        ax.add_patch(rectangle)

    ax.set(xlim=(0, len(mismatch_line)), ylim=(-5, 5))
    fig.suptitle(f"Mutation density distribution of {isolate} {population}, with convolution window of {k}")
    ax.set_xlabel("bp")

    legend_elements = [mpatches.Patch(color="gray", label="no evidence")]
    c = 0
    for name in distributions.keys():  # again cycling in the same order
        legend_elements.append(mpatches.Patch(color="C" + str(c), label=name))
        c += 1

    plt.legend(handles=legend_elements)
    plt.xlim(0, x_axis)
    fig.savefig(out_folder, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--clone", help="path of the clone genome")
    parser.add_argument("--mismatches", help="path of the mismatch array")
    parser.add_argument("--x_axis", help="length of the x axis")
    parser.add_argument("--k", help="convolution window")
    parser.add_argument("--mismatch_threshold", help="threshold for the normalized mismatch density")
    parser.add_argument("--out", help="output path of the plot")

    args = parser.parse_args()
    clone_genome_path = args.clone
    mismatch_array_path = args.mismatches
    x_axis = int(args.x_axis)
    k = int(args.k)
    mismatch_threshold = float(args.mismatch_threshold)
    output_path = args.out

    # set some estetic parameters for the plots
    phage_colors = {"EM11": "C0", "EM60": "C1"}

    population = mismatch_array_path.split("/")[-1].split(".")[0].split("_")[0]
    isolate = mismatch_array_path.split("/")[-1].split(".")[0].split("_")[1]

    len_seq = get_len_seq(clone_genome_path)

    # get all the mismatch arrays resulting from the alignment of the ancestral sequences on the clone genome
    ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)

    # summarise the mismatches of the same ancestral sequence in the same array (normalised). save them to a dictionary.
    mismatch_distributions = summarise_normalise_mismatch_arrays(
        len_seq, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, k
    )

    # create a single array that depicts the normalised mismatch information
    single_array_normalised_mismatches = np.zeros(len_seq)  # this is the mismatch line that we will plot
    c = 0
    for name, array in mismatch_distributions.items():
        # summarise the mismatch density of all phages
        c += 1  # zero for no evidence, 1 for EM11, 2 for EM60 #this order is determined by the order of the lines in the dictionary
        for i in range(len(array)):
            if array[i] >= mismatch_threshold:
                single_array_normalised_mismatches[i] = c

    plot_mismatch_line(
        mismatch_distributions,
        single_array_normalised_mismatches,
        population,
        isolate,
        k,
        x_axis,
        output_path,
    )
