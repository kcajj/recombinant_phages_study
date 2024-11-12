import numpy as np
import csv
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import pysam
from array_compression import decompress_array, retrive_compressed_array_from_str, compress_array
from coordinate_converter import create_coordinate_conversion_map


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


def get_ancestral_alignment(bam_file):

    ancestral_names = []
    ancestral_alignments = []
    ancestral_starts = []
    ancestral_ends = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for ancestral in bam.fetch():
            if not (ancestral.is_secondary):

                ancestral_alignment = []

                mapping_start = ancestral.reference_start
                mapping_end = ancestral.reference_end

                alignment_array = ancestral.get_aligned_pairs()

                for ancestral_pos, clone_pos in alignment_array:
                    if clone_pos != None:  # clips and gaps in the reference are not considered
                        if ancestral_pos == None:
                            ancestral_alignment.append(0)  # gap in the read
                        else:
                            ancestral_alignment.append(1)

                ancestral_names.append(ancestral.query_name)
                ancestral_alignments.append(ancestral_alignment)
                ancestral_starts.append(mapping_start)
                ancestral_ends.append(mapping_end)

    return ancestral_names, ancestral_alignments, ancestral_starts, ancestral_ends


def summarise_ancestral_alignments(
    len_seq,
    ancestral_names,
    ancestral_alignments,
    ancestral_starts,
    ancestral_ends,
):
    """
    create an array that summarises the mapping of ancestral on hybrid
    """

    # create a dictionary with the ancestral names as keys (just one for each ancestral sequence) and the empy mismatch arrays as values,
    distributions = defaultdict(None)
    for i in range(len(ancestral_names)):
        distributions[ancestral_names[i]] = np.zeros(len_seq)

    # extract the mismatches of each ancestral sequence and normalise them
    for name, array in distributions.items():
        for i in range(len(ancestral_names)):
            if ancestral_names[i] == name:
                array = fill(array, ancestral_starts[i], ancestral_ends[i], ancestral_alignments[i])

        distributions[name] = array  # the distributions dictionary stores the normalised arrays

    return distributions


def plot_ancestral_line(
    ancestral_line,
    ancestral_name,
    color,
    ax,
    offset,
    interline,
    thickness,
):

    # compress the array so that it is easy to create rectangles
    plotting_array = compress_array(ancestral_line)

    # plot the rectangles
    height = thickness
    y = offset - interline - thickness
    c = 0

    ax.text(0, y + thickness, f"{ancestral_name}", fontsize=15)

    for lenn, typee in plotting_array:

        x = c
        width = lenn
        c += lenn

        if typee == 0:
            continue  # don't plot gaps

        rectangle = mpatches.Rectangle((x, y), width, height, color=color)
        ax.add_patch(rectangle)

    return y


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


def summarise_normalise_mismatch_arrays(
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
                array = fill(array, mapping_starts[i], mapping_ends[i], mismatch_arrays[i])
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
    mismatch_line,
    population,
    isolate,
    ax,
    offset,
    interline,
    thickness,
    colors,
):
    # compress the array so that it is easy to create rectangles
    plotting_array = compress_array(mismatch_line)

    # plot the rectangles
    height = thickness
    y = offset - interline - thickness
    x_progress = 0

    ax.text(0, y + thickness, f"{population} {isolate}", fontsize=15)

    for lenn, typee in plotting_array:

        x = x_progress
        width = lenn
        x_progress += lenn

        if typee == -1:
            continue  # don't plot gaps

        # invert colors (we are using mimatch density, which is higher in the opposite phage)
        if typee == 1:
            typee = 2
        elif typee == 2:
            typee = 1

        rectangle = mpatches.Rectangle((x, y), width, height, color=colors[typee])
        ax.add_patch(rectangle)

    return y


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--populations", help="population names")
    parser.add_argument("--clones", help="clone names")
    parser.add_argument("--hybrid_ref", help="hybrid reference")
    parser.add_argument("--ancestral_alignment", help="alignment of ancestral sequences on hybridref")
    parser.add_argument("--k", help="convolution window")
    parser.add_argument("--interline", help="interline")
    parser.add_argument("--thickness", help="thickness")
    parser.add_argument("--colors", help="colors")
    parser.add_argument("--out", help="output path of the plot")

    args = parser.parse_args()
    populations = args.populations.split(",")[:-1]
    clones = args.clones.split(",")[:-1]
    hybrid_ref_path = args.hybrid_ref
    ancestral_alignment_path = args.ancestral_alignment
    k = int(args.k)
    interline = int(args.interline)
    thickness = int(args.thickness)
    colors = args.colors.split(",")
    output_path = args.out

    hyb_len = get_len_seq(hybrid_ref_path)

    # get alignments of ancestral sequences.
    ancestral_names, ancestral_alignments, ancestral_starts, ancestral_ends = get_ancestral_alignment(
        ancestral_alignment_path
    )

    ancestral_arrays = summarise_ancestral_alignments(
        hyb_len,
        ancestral_names,
        ancestral_alignments,
        ancestral_starts,
        ancestral_ends,
    )

    n_plots = len(populations) * len(clones) + len(ancestral_arrays.keys())
    fig, ax = plt.subplots(figsize=(20, 10))

    offset = n_plots + (n_plots * (interline + thickness - 1)) + (interline * len(populations))
    y_lim = offset

    c = 1  # take the colors from the second element
    for ancestral_name, array in ancestral_arrays.items():
        offset = plot_ancestral_line(
            array,
            ancestral_name,
            colors[c],
            ax,
            offset,
            interline,
            thickness,
        )
        c += 1

    for population in populations:
        offset -= interline
        for clone in clones:

            mismatch_array_path = f"results/mismatch_arrays/{population}/{population}_{clone}.tsv"
            # recombination_distribution_path=f'results/genomewide_recombination/{population}/{clone}.npz'
            clone_to_hybrid_alignment = f"results/alignments/evidences/{population}/{clone}/{population}_{clone}.bam"
            clone_genome_path = f"data/clones_genomes/{population}/{population}_{clone}.fasta"

            clone_len = get_len_seq(clone_genome_path)

            ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)
            # recombination_distribution = npz_extract(recombination_distribution_path)

            mismatch_distributions = summarise_normalise_mismatch_arrays(
                clone_len, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, k
            )

            # create a single array that simplifies the normalised mismatch information (0 if the normalization is less then one, 1 or 2 in the other case depending on the ancestral reference)
            single_array_normalised_mismatches = np.zeros(clone_len)
            c = 0
            for name, array in mismatch_distributions.items():
                # summarise the mismatch density of all phages
                c += 1  # zero for no evidence, 1 for EM11, 2 for EM60 #this order is determined by the order of the lines in the dictionary!!!!!!!!!! which corresponds tot he order of references in the input file
                for i in range(len(array)):
                    if array[i] == 1:
                        single_array_normalised_mismatches[i] = c

            # convert to the hybrid reference coordinates
            conv_map = create_coordinate_conversion_map(clone_to_hybrid_alignment)
            converted_mismatch_line = np.zeros(hyb_len)
            for clone_coord in range(len(single_array_normalised_mismatches)):
                hyb_coord = conv_map[clone_coord]
                if hyb_coord == None:
                    continue
                converted_mismatch_line[hyb_coord] = single_array_normalised_mismatches[clone_coord]

            # mask out the gaps of the clone on the hybrid reference
            for gap_pos in conv_map["gaps"]:
                converted_mismatch_line[gap_pos] = -1  # gap has value -1 in the mismatch line

            # plot function
            offset = plot_mismatch_line(
                converted_mismatch_line,
                population,
                clone,
                ax,
                offset,
                interline,
                thickness,
                colors,
            )

    # legend
    ax.axes.get_yaxis().set_visible(False)
    ax.set(xlim=(0, hyb_len), ylim=(0, y_lim))
    fig.suptitle(f"Clones. (convolution window {k})", fontsize=20, fontweight="bold")
    plt.savefig(output_path, bbox_inches="tight")
