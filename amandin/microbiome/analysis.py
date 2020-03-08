import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering

# Import Data
abundance_table_path = "./amandin/data/microbiome/zain/abundance_table_100.shared"

with open(abundance_table_path, "r") as file_literal:
    raw_abundance_data = [line.strip().split("\t") for line in file_literal]
    otu_names = raw_abundance_data[0][3:]
    sample_names = [line[1] for line in raw_abundance_data[1:]]
    otu_counts = [line[3:] for line in raw_abundance_data[1:]]

abundance_table = pd.DataFrame(
    data=otu_counts, index=sample_names, columns=otu_names)

# Spectral Clustering
weights = pd.DataFrame(
    squareform(pdist(abundance_table, metric="minkowski", p=1)),
    columns=abundance_table.index,
    index=abundance_table.index)
