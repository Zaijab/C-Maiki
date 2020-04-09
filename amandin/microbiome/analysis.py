# TODO - Scatter Plot Dots by Shape and Color
# TODO - Import GPS Location
# TODO - Scatter plot numbers by site number

# Import Dependencies
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## Instantiation Zone

# File Paths
data_path = "amandin/microbiome/data"
figure_path = "amandin/microbiome/plots"
abundance_table_path = f"{data_path}/abundance_table_100.shared"
metadata_path = f"{data_path}/SuperTransect_mapping_file.csv"

# Abundance Table
with open(abundance_table_path, "r") as file_literal:
    raw_abundance_data = [line.strip().split("\t") for line in file_literal]
    otu_names = raw_abundance_data[0][3:]
    sample_names = list(map(int, [line[1] for line in raw_abundance_data[1:]]))
    otu_counts = [line[3:] for line in raw_abundance_data[1:]]
abundance_table = pd.DataFrame(
    np.array(otu_counts, dtype=np.int64),
    index=sample_names,
    columns=otu_names)
abundance_table["Abundance"] = abundance_table.sum(axis=1)
abundance_table["Presence"] = abundance_table.drop("Abundance", axis=1).where(
    abundance_table == 0, 1).sum(axis=1)

# Metadata
metadata = pd.read_csv(metadata_path, index_col = 0)

# Only look at Plantroot groups
plantroot_ids = metadata.loc[metadata["collection_label"] == "PlantRoot"]
plantroot_group = abundance_table.filter(items = list(plantroot_ids.index), axis = 0)
groups = [plantroot_group]

# Compute Distance Matrix
group_weight = [
    pd.DataFrame(
        squareform(pdist(group.drop(["Abundance", "Presence"], axis=1),
                         metric="minkowski", p=1)),
        columns=group.index,
        index=group.index)
    for group in groups
]

# Normalize Distance Matrix
normalized_group_weight = [
    weight / np.amax(weight) for weight in group_weight]

# Compute Kernel exp(-(w_ij)^2/sigma^2)
group_kernel = [np.exp(- (weight - np.mean(weight)) ** 2 / np.var(weight))
                for weight in normalized_group_weight]

# Compute Diagonal Matrix
group_diagonal = [
    pd.DataFrame(
        np.diag(np.sum(kernel)),
        columns=kernel.columns,
        index=kernel.index)
    for kernel in group_kernel]

# Compute Laplacian
group_laplacian = [
    pd.DataFrame(
        diagonal - weight,
        columns=diagonal.columns,
        index=diagonal.index)
    for diagonal, weight in zip(group_diagonal, group_kernel)]

# Compute Eigenvectors / Eigenvalues
group_eigen = [np.linalg.eig(laplacian)
               for laplacian in group_laplacian]

# Plot
for i, eigen in enumerate(group_eigen):
    eigenvalues, eigenvectors = eigen
    eigenvectors = eigenvectors.real
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(eigenvectors[1], eigenvectors[2], eigenvectors[3])
    plt.savefig(f"{figure_path}/two_eigenvectors_non_normalized_plantroot3d.pdf")
    plt.close()
