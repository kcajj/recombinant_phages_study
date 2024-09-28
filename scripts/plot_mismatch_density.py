import numpy as np
import csv
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict
from array_compression import decompress_array, retrive_compressed_array_from_str

def get_mismatch_arrays(mismatch_array_path):
    csv.field_size_limit(sys.maxsize)

    ancestral_names=[]
    mismatch_arrays=[]
    mapping_starts=[]
    mapping_ends=[]

    with open(mismatch_array_path) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            #extracting mismatch arrays and mapping info
            ancestral_name=line[0]
            mapping_start=int(line[1])
            mapping_end=int(line[2])

            compressed_mismatch_array=retrive_compressed_array_from_str(line[3])
            mismatch_array=decompress_array(compressed_mismatch_array)

            ancestral_names.append(ancestral_name)
            mismatch_arrays.append(mismatch_array)
            mapping_starts.append(mapping_start)
            mapping_ends.append(mapping_end)

    return ancestral_names, mismatch_arrays, mapping_starts, mapping_ends

def npz_extract(npz_file):
    npz=np.load(npz_file)
    lst=npz.files
    for item in lst:
        array=npz[item]
    return array

def get_len_seq(fasta_path):
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return int(len(record.seq))

def plot_mismatches(clone_genome_path, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, population, isolate, k, out_folder):
    fig=plt.figure(figsize=(18,6))
    
    len_seq = get_len_seq(clone_genome_path)
    
    lines=defaultdict(None)
    for i in range(len(ancestral_names)):
        lines[ancestral_names[i]]=np.zeros(len_seq)

    for name,array in lines.items():
        for i in range(len(ancestral_names)):
            if ancestral_names[i]==name:
                for j in range(len(array)):
                    if j>=mapping_starts[i] and j<=mapping_ends[i]:
                        array[j]+=mismatch_arrays[i][j-mapping_starts[i]]
        convolved_mismatch_array = np.convolve(array, np.ones(k), 'valid') / k
        plt.plot(convolved_mismatch_array,label=name,alpha=0.5,color=phage_colors[name])
    
    plt.legend()
    plt.title(f'Mutation density distribution of {isolate} {population}, with convolution window of {k}')
    plt.ylabel('mutation density')
    plt.xlabel('bp')
    fig.savefig(out_folder, bbox_inches='tight')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--clone", help="path of the clone genome")
    parser.add_argument("--mismatches", help="path of the mismatch array")
    parser.add_argument("--out", help="output path of the plot")

    args = parser.parse_args()
    clone_genome_path=args.clone
    mismatch_array_path=args.mismatches
    output_path=args.out

    k=500 #convolution window

    #set some estetic parameters for the plots
    phage_colors={'EM11':'C0','EM60':'C1','EC2D2':'C2'}

    population=mismatch_array_path.split("/")[-1].split(".")[0].split("_")[0]
    isolate=mismatch_array_path.split("/")[-1].split(".")[0].split("_")[1]
    #plot the clone assemblies mutation density
    
    ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)
    
    plot_mismatches(clone_genome_path, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, population, isolate, k, output_path)