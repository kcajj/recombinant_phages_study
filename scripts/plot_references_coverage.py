import numpy as np
from handle_msa import length_msa, extract_references_names
import csv
import sys
import matplotlib.pyplot as plt
import plotly.express as px
from hmm_prediction_arrays import get_evidence_arrays
from array_compression import npz_extract


def plot_references_coverage(array1, array2, output_path, title, x_label, y_label, refs_msa_path, population):

    ### temporary: adding the evidence arrays to the plot
    evidences_file = f"results/evidence_arrays/{population}/{id}.tsv"
    read_names, evidence_arrays, mapping_starts, mapping_ends, c_reads = get_evidence_arrays(evidences_file)
    evidences_to_plot = []
    x = []
    for i in range(len(read_names)):
        if read_names[i] == "P1_C1":
            for pos, val in enumerate(evidence_arrays[i]):
                if val != 0:
                    x.append(pos + mapping_starts[i])
                if val == 1:
                    evidences_to_plot.append("blue")
                elif val == 2:
                    evidences_to_plot.append("orange")
    y = [1] * len(x)
    ###

    plt.figure(figsize=(20, 5))
    plt.scatter(x, y, marker="|", c=evidences_to_plot, alpha=0.3)
    plt.plot(array1, alpha=0.5)
    plt.plot(array2, alpha=0.5)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(["evidences"] + extract_references_names(refs_msa_path))
    plt.savefig(output_path)

    html_path = output_path[:-3] + "html"
    fig = px.line(x=range(len(array1)), y=array1)
    fig.add_scatter(x=np.array(range(len(array2))), y=array2, mode="lines")
    fig.write_html(html_path)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Makes a prediction on each evidence array and summarises the information in a single recombination array",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--msa_refs", help="path of the hybrid reference")
    parser.add_argument("--coverage", help="path of the folder containing coverage arrays")
    parser.add_argument("--out", help="output path of the .npz file containing the coverage plot")

    args = parser.parse_args()
    refs_msa_path = args.msa_refs
    coverage_path = args.coverage
    output_path = args.out

    coverage_0_path = f"{coverage_path[:-4]}_0.npz"
    coverage_1_path = f"{coverage_path[:-4]}_1.npz"

    coverage = npz_extract(coverage_path)
    coverage0 = npz_extract(coverage_0_path)
    coverage1 = npz_extract(coverage_1_path)

    id = output_path.split("/")[-1].split(".")[0]
    population = output_path.split("/")[-1].split(".")[0].split("_")[0]

    plot_references_coverage(
        coverage0,
        coverage1,
        output_path,
        f"References coverage {id}",
        "position",
        "coverage",
        refs_msa_path,
        population,
    )

    normalised0 = np.divide(
        coverage0.astype(float),
        coverage.astype(float),
        out=np.zeros_like(coverage0.astype(float)),
        where=coverage.astype(float) != 0,
    )
    normalised1 = np.divide(
        coverage1.astype(float),
        coverage.astype(float),
        out=np.zeros_like(coverage1.astype(float)),
        where=coverage.astype(float) != 0,
    )

    normalised_output_path = output_path[:-4] + "_normalised.pdf"

    plot_references_coverage(
        normalised0,
        normalised1,
        normalised_output_path,
        f"Normalised references coverage {id}",
        "position",
        "coverage rate",
        refs_msa_path,
        population,
    )
