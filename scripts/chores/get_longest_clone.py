from Bio import SeqIO

if __name__ == "__main__":
    populations = ["P1", "P2", "P3"]
    clones = ["C1", "C2", "C3", "C4"]

    max_length = 0
    for population in populations:
        for clone in clones:
            path = f"data/clones_genomes/{population}/{population}_{clone}.fasta"
            for record in SeqIO.parse(path, "fasta"):
                length = len(record.seq)
                if length > max_length:
                    max_length = length
                    longest_clone = f"{population}_{clone}"

    print(longest_clone)
    print(max_length)
