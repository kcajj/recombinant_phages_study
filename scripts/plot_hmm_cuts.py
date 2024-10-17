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


def plot_hmm_cuts(
    clone_genome_path,
    ancestral_names,
    mismatch_arrays,
    mapping_starts,
    mapping_ends,
    recombination_distribution,
    population,
    isolate,
    k,
    x_axis,
    out_folder,
):
    fig = plt.figure(figsize=(18, 6))

    len_seq = get_len_seq(clone_genome_path)

    lines = defaultdict(None)
    for i in range(len(ancestral_names)):
        lines[ancestral_names[i]] = np.zeros(len_seq)

    total_mismatches = np.zeros(len_seq)
    for i in range(len(ancestral_names)):
        for j in range(len(mismatch_arrays[i])):
            if j >= mapping_starts[i] and j <= mapping_ends[i]:
                total_mismatches[j] += mismatch_arrays[i][j - mapping_starts[i]]

    both_mismatches = np.zeros(len_seq)

    c = 0  # zero for no evidence, 1 for EM11, 2 for EM60
    for name, array in lines.items():
        for i in range(len(ancestral_names)):
            if ancestral_names[i] == name:
                for j in range(len(array)):
                    if j >= mapping_starts[i] and j <= mapping_ends[i]:
                        array[j] += mismatch_arrays[i][j - mapping_starts[i]]
        convolved_mismatch_array = np.convolve(array, np.ones(k), "valid") / k
        convolved_total_array = np.convolve(total_mismatches, np.ones(k), "valid") / k
        normalised_array = np.divide(
            convolved_mismatch_array,
            convolved_total_array,
            out=np.zeros_like(convolved_mismatch_array),
            where=convolved_total_array != 0,
        )

        # summarise the mismatch density of all phages
        c += 1  # zero for no evidence, 1 for EM11, 2 for EM60 #this order is determined by the order of the lines in the dictionary
        for i in range(len(ancestral_names)):
            if ancestral_names[i] == name:
                for j in range(len(normalised_array)):
                    if normalised_array[j] == 1:
                        both_mismatches[j] = c

    # create a colour array with 0=gray and from 1 onwards C0, C1, C2, C3 ecc
    colors = {"C0": "gray"}
    c = 0
    for line in lines.keys():
        colors["C" + str(c + 1)] = "C" + str(c)
        c += 1
        # cycling through the lines in the same order as we did in the previous loop, to guarantee correct colours (if a line is first it gets c0)

    fig, ax = plt.subplots(nrows=1, figsize=(18, 6))
    plotting_array = compress_array(both_mismatches)

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
    ax.set(xlim=(0, len(both_mismatches)), ylim=(-5, 5))

    # change coordinates
    converted_recombination = np.zeros_like(recombination_distribution)
    for hybrid_coord in range(len(recombination_distribution)):
        clone_coord = coordinate_conversion_map[hybrid_coord]
        if clone_coord == None:
            continue
        converted_recombination[clone_coord] = recombination_distribution[hybrid_coord]
    cuts = [i for i, value in enumerate(converted_recombination) if value == 1]
    ax.scatter(cuts, [0] * len(cuts), color="red", marker="|", s=1200)

    fig.suptitle(f"Mutation density distribution of {isolate} {population}, with convolution window of {k}")
    ax.set_xlabel("bp")

    legend_elements = [mpatches.Patch(color="gray", label="no evidence")]
    c = 0
    for line in lines.keys():  # again cycling in the same order
        legend_elements.append(mpatches.Patch(color="C" + str(c), label=line))
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
    parser.add_argument("--recombination", help="path of the .npz recombination array")
    parser.add_argument("--hc_alignment", help="path of the hybrid-clone alignment")
    parser.add_argument("--x_axis", help="length of the x axis")
    parser.add_argument("--k", help="convolution window")
    parser.add_argument("--out", help="output path of the plot")

    args = parser.parse_args()
    clone_genome_path = args.clone
    mismatch_array_path = args.mismatches
    recombination_distribution_path = args.recombination
    hybrid_clone_alignment = args.hc_alignment
    x_axis = int(args.x_axis)
    k = int(args.k)
    output_path = args.out

    # set some estetic parameters for the plots
    phage_colors = {"EM11": "C0", "EM60": "C1"}

    population = mismatch_array_path.split("/")[-1].split(".")[0].split("_")[0]
    isolate = mismatch_array_path.split("/")[-1].split(".")[0].split("_")[1]
    # plot the clone assemblies mutation density

    ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)

    recombination_distribution = npz_extract(recombination_distribution_path)

    coordinate_conversion_map = get_hyb_ref_map(hybrid_clone_alignment)

    plot_hmm_cuts(
        clone_genome_path,
        ancestral_names,
        mismatch_arrays,
        mapping_starts,
        mapping_ends,
        recombination_distribution,
        population,
        isolate,
        k,
        x_axis,
        output_path,
    )
