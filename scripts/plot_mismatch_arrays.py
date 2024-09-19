import numpy as np
import csv
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import defaultdict
from array_compression import decompress_array

def get_mismatch_arrays(mismatch_array_path):
    csv.field_size_limit(sys.maxsize)

    ancestral_names=[]
    mismatch_arrays=[]
    mapping_starts=[]
    mapping_ends=[]

    with open(mismatch_array_path) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            ancestral_name=line[0]
            mapping_start=int(line[1])
            mapping_end=int(line[2])

            mismatch_array_str_tuples=line[3][2:-2].split('), (')

            compressed_mismatch_array = []
            for tuple_str in mismatch_array_str_tuples:
                x,y=tuple_str.split(', ')
                compressed_mismatch_array.append((int(x),int(y)))
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

'''
def plot_mappings(mismatch_array, coverage_array, population, isolate, k, out_folder):
    
    #plot the distribution of mismatches of a bam file
    
    number_of_plots=len(mismatch_distribution.keys())+1
    fig, axs =plt.subplots(number_of_plots,sharex=True,constrained_layout = True, figsize=(8,10))

    fig.suptitle(f'mutation density distribution between {population}-{isolate} and references, with convolution window of {k}')

    for reference, distribution in mismatch_distribution.items():
        distribution=np.convolve(distribution,np.ones(k),'valid')/k
        l=len(distribution)
        x=np.linspace(0,l,l)
        axs[0].plot(x,distribution,color=phage_colors[reference])
    axs[0].legend(mismatch_distribution.keys())
    axs[0].set_title(reference)
    axs[0].set_ylabel('mutation density')
    axs[0].set_xlabel('bp')

    #plot coverage
        
    fig.savefig(out_folder, bbox_inches='tight')
    plt.close(fig)
'''

def plot_mismatches(clone_genome_path, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, coverage_array, population, isolate, k, out_folder):
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
        array = np.convolve(array, np.ones(k), 'valid') / k
        plt.plot(array,label=name,alpha=0.5)
    
    plt.legend()
    plt.title(f'Mutation density distribution of {isolate} {population}, with convolution window of {k}')
    plt.ylabel('mutation density')
    plt.xlabel('bp')
    fig.savefig(out_folder, bbox_inches='tight')

if __name__ == "__main__":
    
    populations=['P1']
    isolates=['C1']
    k=1000 #convolution window

    #set some estetic parameters for the plots
    phage_colors={'EC2D2':'C2','EM11':'C0','EM60':'C1'}
    mapping_colors={'primary':'g','supplementary':'b','secondary':'r'}

    #plot the clone assemblies mutation density
    for population in populations:
        for isolate in isolates:
            
            clone_genome_path=f'data/clones_genomes/{population}/{population}_{isolate}.fasta'
            mismatch_array_path = f'results/clones/mismatch_arrays/{population}/{population}_{isolate}.tsv'
            coverage_array_path = f'results/clones/coverage_arrays/{population}/{population}_{isolate}.npz'
            out_folder = f'results/clones/plots/{population}_{isolate}.png'

            ancestral_names, mismatch_arrays, mapping_starts, mapping_ends = get_mismatch_arrays(mismatch_array_path)
            coverage_array=decompress_array(npz_extract(coverage_array_path))
            
            plot_mismatches(clone_genome_path, ancestral_names, mismatch_arrays, mapping_starts, mapping_ends, coverage_array, population, isolate, k, out_folder)