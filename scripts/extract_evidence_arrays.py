import numpy as np
from handle_msa import length_msa, get_evidences_distributions, add_to_msa
import pysam
import time
import csv
from array_compression import compress_array
import matplotlib.pyplot as plt


def write_evidence_arrays(bam_file, refs_msa_path, output_evidences_path):

    time_spent_per_read = []
    tot_time_start = time.time()
    c_tot_alignments = 0
    c_useful_alignments = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_evidences_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
            for read in bam.fetch():
                if not (read.is_secondary):

                    start_time = time.time()

                    read_seq = read.query_sequence
                    read_msa_seq = ""

                    mapping_start = read.reference_start
                    mapping_end = read.reference_end

                    alignment_array = read.get_aligned_pairs()

                    for read_pos, ref_pos in alignment_array:
                        if ref_pos != None:  # clips and gaps in the reference are not considered
                            if read_pos == None:
                                read_msa_seq += "-"  # gap in the read
                            else:
                                read_msa_seq += read_seq[read_pos].lower()

                    msa_matrix = add_to_msa(refs_msa_path, read_msa_seq, mapping_start, mapping_end)

                    e_distribution_to_plot = get_evidences_distributions(msa_matrix)

                    e_distribution = np.where(
                        e_distribution_to_plot > 0, e_distribution_to_plot - 1, e_distribution_to_plot
                    )
                    compressed_e_distribution = compress_array(e_distribution)

                    np.set_printoptions(threshold=np.inf, linewidth=np.inf)
                    writer.writerow([read.query_name, mapping_start, mapping_end, compressed_e_distribution])

                    end_time = time.time()
                    time_spent_per_read.append(end_time - start_time)

                    c_useful_alignments += 1
                    """
                    plot_path = f"results/ev_plots/{read.query_name}_{int(time.time()*1000)}.png"

                    hmm_plot, (evidences, prediction) = plt.subplots(2, 1, figsize=(10, 5))
                    hmm_plot.suptitle(f'HMM read {c_tot_alignments}')

                    colours = np.where(e_distribution_to_plot == 0, "green", np.where(e_distribution_to_plot == 1, "red", np.where(e_distribution_to_plot == 2, "blue", "orange")))
                    evidences.scatter(range(mapping_start,len(e_distribution_to_plot)+mapping_start), e_distribution_to_plot, c=colours, marker='|', alpha=0.5)
                    evidences.set_title(f'evidence distribution (0=same, 1=error, 2=evidence for)')
                    evidences.set_xlabel("basepair")
                    evidences.set_ylabel("visible states")
                    hmm_prediction = np.zeros(len(e_distribution_to_plot))

                    colours = np.where(hmm_prediction == 0, "blue", "orange")
                    prediction.scatter(range(mapping_start,len(e_distribution_to_plot)+mapping_start), hmm_prediction, c=colours, marker='|', alpha=0.5)
                    prediction.set_title(f'HMM prediction)')
                    prediction.set_xlabel("basepair")
                    prediction.set_ylabel("hidden states")

                    hmm_plot.tight_layout()
                    hmm_plot.savefig(plot_path)
                    plt.close(hmm_plot)
                    """

                c_tot_alignments += 1

    output_stats_path = output_evidences_path[:-4] + "_stats.txt"
    with open(output_stats_path, "w") as f:
        f.write("evidence arrays extraction run of " + bam_file + "\n")
        f.write("mean time spent per read " + str(np.mean(time_spent_per_read)) + "\n")
        f.write("total time " + str(time.time() - tot_time_start) + "\n")
        f.write("total reads " + str(c_tot_alignments) + "\n")
        f.write("reads used " + str(c_useful_alignments) + "\n")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Extract all the evidence arrays from the alignments of a bam file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", help="path of the bam file")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--evidences_out", help="output path of the .tsv file containing the evidence arrays")

    args = parser.parse_args()
    bam_file = args.bam
    refs_msa_path = args.msa_refs
    output_evidences_path = args.evidences_out

    write_evidence_arrays(bam_file, refs_msa_path, output_evidences_path)
