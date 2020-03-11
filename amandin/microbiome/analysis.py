# Import Data
#################################################################
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import SpectralClustering
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_path = "amandin/microbiome/data"
figure_path = "amandin/microbiome/plots"
abundance_table_path = f"{data_path}/abundance_table_100.shared"

with open(abundance_table_path, "r") as file_literal:
    raw_abundance_data = [line.strip().split("\t") for line in file_literal]
    otu_names = raw_abundance_data[0][3:]
    sample_names = [line[1] for line in raw_abundance_data[1:]]
    otu_counts = [line[3:] for line in raw_abundance_data[1:]]

abundance_table = pd.DataFrame(
    np.array(otu_counts, dtype=np.int64),
    index=sample_names,
    columns=otu_names)

abundance_table["Abundance"] = abundance_table.sum(axis=1)
abundance_table["Presence"] = abundance_table.drop(
    'Abundance', axis=1).where(abundance_table > 0, 1).sum(axis=1)

#################################################################

print(abundance_table)
