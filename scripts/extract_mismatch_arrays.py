import numpy as np
import pysam
import csv
from Bio import SeqIO
from array_compression import compress_array


def get_seq(fasta_path):
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return str(record.seq)


def write_mismatches_arrays(clone_genome_path, bam_file, output_mismatch_path):

    clone_seq = get_seq(clone_genome_path)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_mismatch_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
            for ancestral in bam.fetch():
                if not (ancestral.is_secondary):

                    ancestral_seq = ancestral.query_sequence
                    ancestral_full_seq = ""
                    mismatches = np.zeros(len(clone_seq), dtype=int)

                    mapping_start = ancestral.reference_start
                    mapping_end = ancestral.reference_end

                    alignment_array = ancestral.get_aligned_pairs()

                    for ancestral_pos, clone_pos in alignment_array:
                        if clone_pos != None:  # clips and gaps in the reference are not considered
                            if ancestral_pos == None:
                                ancestral_full_seq += "-"  # gap in the read
                            else:
                                ancestral_full_seq += ancestral_seq[ancestral_pos].upper()

                    clone_seq_crop = clone_seq[mapping_start:mapping_end]

                    for i in range(len(clone_seq_crop)):
                        if ancestral_full_seq[i] != "-" and clone_seq_crop[i] != "-":
                            if ancestral_full_seq[i] != clone_seq_crop[i]:
                                mismatches[i] += 1

                    compressed_mismatches = compress_array(mismatches)

                    np.set_printoptions(threshold=np.inf, linewidth=np.inf)
                    writer.writerow([ancestral.query_name, mapping_start, mapping_end, compressed_mismatches])


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Extract the mismatch arrays from the alignments of a bam file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", help="path of the bam file")
    parser.add_argument("--clone", help="path of the clone genome")
    parser.add_argument("--mismatch_out", help="output path of the .tsv file containing the evidence arrays")

    args = parser.parse_args()
    bam_file = args.bam

    clone_genome_path = args.clone
    output_mismatch_path = args.mismatch_out

    write_mismatches_arrays(clone_genome_path, bam_file, output_mismatch_path)
