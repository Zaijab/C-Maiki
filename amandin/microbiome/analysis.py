# TODO - Scatter Plot Dots by Shape and Color
# TODO - Scatter Plot in 3D with 3 eigenvectors
# TODO - Import GPS Location
# TODO - Scatter plot numbers by site number

from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_path = "amandin/microbiome/data"
figure_path = "amandin/microbiome/plots"
abundance_table_path = f"{data_path}/abundance_table_100.shared"
metadata_path = f"{data_path}/SuperTransect_mapping_file.csv"

with open(abundance_table_path, "r") as file_literal:
    raw_abundance_data = [line.strip().split("\t") for line in file_literal]
    otu_names = raw_abundance_data[0][3:]
    sample_names = [line[1] for line in raw_abundance_data[1:]]
    otu_counts = [line[3:] for line in raw_abundance_data[1:]]

with open(metadata_path, "r") as file_literal:
    raw_metadata = [line.strip().split(",") for line in file_literal][1:]
metadata = {}
for line in raw_metadata:
    metadata[line[0]] = line[20]


abundance_table = pd.DataFrame(
    np.array(otu_counts, dtype=np.int64),
    index=sample_names,
    columns=otu_names)

abundance_table["Abundance"] = abundance_table.sum(axis=1)
abundance_table["Presence"] = abundance_table.drop("Abundance", axis=1).where(
    abundance_table == 0, 1).sum(axis=1)


plantroot_ids = [key for key, value in metadata.items() if value == "PlantRoot"] 
plantroot_group = abundance_table.filter(items = plantroot_ids, axis = 0)

groups = [plantroot_group]

group_weight = [
    pd.DataFrame(
        squareform(pdist(group.drop(["Abundance", "Presence"], axis=1),
                         metric="minkowski", p=1)),
        columns=group.index,
        index=group.index)
    for group in groups
]

normalized_group_weight = [
    weight / np.amax(weight) for weight in group_weight]

group_kernel = [np.exp(- (weight - np.mean(weight)) ** 2 / np.var(weight))
                for weight in normalized_group_weight]

group_diagonal = [
    pd.DataFrame(
        np.diag(np.sum(kernel)),
        columns=kernel.columns,
        index=kernel.index)
    for kernel in group_kernel]


group_laplacian = [
    pd.DataFrame(
        diagonal - weight,
        columns=diagonal.columns,
        index=diagonal.index)
    for diagonal, weight in zip(group_diagonal, group_kernel)]

group_weight_hat = [
    pd.DataFrame(
        np.linalg.inv(diagonal) * weight,
        columns=diagonal.columns,
        index=diagonal.index)
    for diagonal, weight in zip(group_diagonal, normalized_group_weight)]


group_eigen = [np.linalg.eig(laplacian)
               for laplacian in group_laplacian]


for i, eigen in enumerate(group_eigen):
    eigenvalues, eigenvectors = eigen
    print(eigenvectors[1],eigenvectors[2],eigenvectors[3])
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    plt.title(f"Plantroot")
    ax.scatter3D(eigenvectors[1], eigenvectors[2], eigenvectors[3])
    plt.savefig(f"{figure_path}/two_eigenvectors_non_normalized_plantroot3d.pdf")
    plt.close()
