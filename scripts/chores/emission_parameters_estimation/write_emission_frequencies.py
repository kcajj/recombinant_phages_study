import numpy as np
import pysam
import time
from Bio import AlignIO


def get_evidences_distributions(msa_matrix, i_ref1=0, i_ref2=1, i_extra=2):
    l = len(msa_matrix[0])
    e_distribution = np.zeros(l, dtype=int)

    for pos, array in enumerate(msa_matrix[0]):
        nuc_extra = msa_matrix[i_extra, pos]
        nuc_first_ref = msa_matrix[i_ref1, pos]
        nuc_second_ref = msa_matrix[i_ref2, pos]
        # if nuc_first_ref!='-' and nuc_second_ref!='-':
        if nuc_extra != "-" and nuc_first_ref != "-" and nuc_second_ref != "-":
            if nuc_extra == nuc_first_ref and nuc_extra == nuc_second_ref:
                continue
            elif nuc_extra != nuc_first_ref and nuc_extra != nuc_second_ref:
                e_distribution[pos] = 1
            elif nuc_extra == nuc_first_ref and nuc_extra != nuc_second_ref:
                e_distribution[pos] = 2
            elif nuc_extra != nuc_first_ref and nuc_extra == nuc_second_ref:
                e_distribution[pos] = 3
        # elif nuc_first_ref=="-" and nuc_second_ref=="-":
        #    if nuc_extra==nuc_first_ref:
        #        continue
        #    elif nuc_extra!=nuc_first_ref:
        #        e_distribution[pos]=1

    return e_distribution


def add_to_msa(msa_path, seq, mapping_start, mapping_end):
    alignment = AlignIO.read(open(msa_path), "fasta")
    l = alignment.get_alignment_length()

    msa_matrix = np.zeros([3, l], dtype=str)
    for i, record in enumerate(alignment):
        for pos, nuc in enumerate(record.seq):
            msa_matrix[i][pos] = nuc

    cut_msa_matrix = msa_matrix[:, mapping_start:mapping_end]
    for pos in range(len(cut_msa_matrix[2])):
        cut_msa_matrix[2][pos] = seq[pos]

    return cut_msa_matrix


def get_emission_frequencies(bam_file, refs_msa_path):

    null_freq = []
    a_freq = []
    b_freq = []

    time_spent_per_read = []
    c_tot_alignments = 0
    c_useful_alignments = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not (read.is_secondary):

                start_time = time.time()

                name = read.query_name

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

                l = len(e_distribution)
                null_freq.append(np.count_nonzero(e_distribution == 0) / l)
                a_freq.append(np.count_nonzero(e_distribution == 1) / l)
                b_freq.append(np.count_nonzero(e_distribution == 2) / l)

                end_time = time.time()
                time_spent_per_read.append(end_time - start_time)

                c_useful_alignments += 1

            c_tot_alignments += 1

    mean_time = np.mean(time_spent_per_read)

    return null_freq, a_freq, b_freq, mean_time, c_tot_alignments, c_useful_alignments


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Extract all the evidence arrays from the alignments of a bam file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", help="path of the bam file of the pure phages mapped on the hybrid reference")
    parser.add_argument("--msa_refs", help="path of the msa between the references")
    parser.add_argument("--name", help="name of the phage")
    parser.add_argument("--out", help="output path of the file containing the emission frequencies")

    args = parser.parse_args()
    bam_file = args.bam
    refs_msa_path = args.msa_refs
    name = args.name
    out = args.out

    tot_time_start = time.time()

    null_freq, a_freq, b_freq, mean_time, c_tot_alignments, c_useful_alignments = get_emission_frequencies(
        bam_file, refs_msa_path
    )

    with open(out, "a") as f:
        f.write(name + "\n")
        f.write("\n")
        f.write("null_freq " + str(np.mean(null_freq)) + "\n")
        f.write("a_freq " + str(np.mean(a_freq)) + "\n")
        f.write("b_freq " + str(np.mean(b_freq)) + "\n")
        f.write("\n")
        f.write("mean time spent per read " + str(mean_time) + "\n")
        f.write("total time " + str(time.time() - tot_time_start) + "\n")
        f.write("total reads " + str(c_tot_alignments) + "\n")
        f.write("reads used " + str(c_useful_alignments) + "\n")
        f.write("\n")
