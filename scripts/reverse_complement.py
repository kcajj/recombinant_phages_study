import sys
from Bio import SeqIO
from Bio.Seq import Seq

#i used this script to reverse complement the clone genomes that were mapping inverse to the reference genomes

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def process_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    
    for record in records:
        record.seq = record.seq.reverse_complement()
    
    SeqIO.write(records, file_path, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python reverse_complement.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    process_fasta(fasta_file)
