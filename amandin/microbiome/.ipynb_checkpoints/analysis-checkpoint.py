# Import Dependencies
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import eigs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Instantiation / Data Cleaning Zone
# File Paths
data_path = "amandin/microbiome/data"
figure_path = "amandin/microbiome/plots"
abundance_table_path = f"{data_path}/abundance_table_97.shared"
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
metadata = pd.read_csv(metadata_path, index_col=0)

mosquito_metadata = metadata.loc[metadata["sample_type"] == "Mosquito"]
mosquito_abundance = abundance_table.filter(
    items=list(mosquito_metadata.index), axis=0)


mosquito_adj = squareform(pdist(mosquito_abundance.drop(
    ["Abundance", "Presence"], axis=1).drop([105525, 105502], axis=0), metric="minkowski", p=1))

kernel = np.exp(- (mosquito_adj ** 2) / (3000**2))

diagonal = np.diag(np.sum(kernel, axis=1))

laplacian = diagonal - mosquito_adj

#normalized_laplacian = np.linalg.inv(diagonal) * laplacian

eigenvalues, eigenvectors = eigs(laplacian, k=4)
eigenvectors = eigenvectors.T.real

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter3D(eigenvectors[1], eigenvectors[2], eigenvectors[3])
plt.title("Generalized Eigenvector Mosquito 3D")
plt.savefig(f"{figure_path}/generalized_eigenvector_mosquito_3d.png")
plt.close()
