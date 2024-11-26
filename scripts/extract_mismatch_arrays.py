import csv

import numpy as np
import pysam
from array_compression import compress_array
from Bio import SeqIO


def get_seq(fasta_path):
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return str(record.seq)


def write_mismatches_arrays(clone_genome_path, bam_file, output_mismatch_path):
    """
    Extracts mismatch information from alignments of ancestral sequences to a clone genome
    and writes the mismatch data to a specified output file.

    Parameters:
    - clone_genome_path: Path to the clone genome fasta file.
    - bam_file: Path to the BAM file containing alignments.
    - output_mismatch_path: Path to the output file where mismatch arrays will be written.
    """

    # Retrieve the clone sequence from the provided fasta file
    clone_seq = get_seq(clone_genome_path)

    # Open the BAM file and output TSV for writing mismatch data
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_mismatch_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")

            # Iterate over each ancestral sequence in the BAM file
            for ancestral in bam.fetch():
                if not (ancestral.is_secondary):

                    # Retrieve the ancestral sequence and initialize variables
                    ancestral_seq = ancestral.query_sequence
                    ancestral_full_seq = ""
                    mismatches = np.zeros(len(clone_seq), dtype=int)

                    # Get the start and end positions of the alignment on the clone
                    mapping_start = ancestral.reference_start
                    mapping_end = ancestral.reference_end

                    # Retrieve aligned pairs of positions between ancestral and clone
                    alignment_array = ancestral.get_aligned_pairs()

                    # Construct the full ancestral sequence considering gaps
                    for ancestral_pos, clone_pos in alignment_array:
                        if clone_pos != None:
                            if ancestral_pos == None:
                                ancestral_full_seq += "-"
                            else:
                                ancestral_full_seq += ancestral_seq[ancestral_pos].upper()

                    # Extract the corresponding segment from the clone sequence
                    clone_seq_crop = clone_seq[mapping_start:mapping_end]

                    # Compare ancestral and clone sequences to identify mismatches
                    for i in range(len(clone_seq_crop)):
                        if clone_seq_crop[i] != "-":
                            if ancestral_full_seq[i] != clone_seq_crop[i]:
                                mismatches[i] += 1
                        else:
                            mismatches[i] += 1

                    # Compress and write the mismatch information to the output file
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
