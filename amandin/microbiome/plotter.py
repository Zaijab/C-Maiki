presence_count = np.sum(
    np.where(abundance_table > 0, 1, 0), axis=1)
abundance_count = np.sum(abundance_table, axis=1)
table = np.column_stack((presence_count, abundance_count))
table_leq_20 = table[table[:, 0] < 20]
table_ged_20_leq_600 = table[np.logical_and(
    20 <= table[:, 0], table[:, 0] < 600)]
table_ged_600 = table[600 <= table[:, 0]]

# Precence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples All")
plt.hist(table[:, 0], bins=range(0, 900, 20))
plt.savefig(f"{figure_path}/presence_histogram_group1.pdf")
plt.close()

# Abundance Count Histogram
plt.figure()
plt.xlabel("Abundance")
plt.ylabel("Samples")
plt.title("Abundance vs Samples All")
plt.hist(table[:, 1])
plt.savefig(f"{figure_path}/abunance_histogram_group1.pdf")
plt.close()

# Abundance vs Presence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance All")
plt.scatter(x=table[:, 0], y=table[:, 1])
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_group1.pdf")
plt.close()

# Less than 20

# Precence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples Less Than 20")
plt.hist(table_leq_20[:, 0], bins=range(0, 20, 1))
plt.savefig(f"{figure_path}/presence_histogram_group2.pdf")
plt.close()

# Abundance Count Histogram
plt.figure()
plt.xlabel("Abundance")
plt.ylabel("Samples")
plt.title("Abundance vs Samples Less Than 20")
plt.hist(table_leq_20[:, 1])
plt.savefig(f"{figure_path}/abunance_histogram_group2.pdf")
plt.close()

# Abundance vs Presence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance Less Than 20")
plt.scatter(x=table_leq_20[:, 0], y=table_leq_20[:, 1])
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_group2.pdf")
plt.close()


# Less than 20

# Precence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples Greater than or Equal to 20 Less than 600")
plt.hist(table_ged_20_leq_600[:, 0], bins=range(20, 610, 10))
plt.savefig(f"{figure_path}/presence_histogram_group3.pdf")
plt.close()

# Abundance Count Histogram
plt.figure()
plt.xlabel("Abundance")
plt.ylabel("Samples")
plt.title("Abundance vs Samples Greater than or Equal to 20 Less than 600")
plt.hist(table_ged_20_leq_600[:, 1])
plt.savefig(f"{figure_path}/abunance_histogram_group3.pdf")
plt.close()

# Abundance vs Presence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance Greater than or Equal to 20 Less than 600")
plt.scatter(x=table_ged_20_leq_600[:, 0], y=table_ged_20_leq_600[:, 1])
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_group3.pdf")
plt.close()

# Less than 20

# Precence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples Greater than or Equal to 600")
plt.hist(table_ged_600[:, 0], bins=range(590, 900, 10))
plt.savefig(f"{figure_path}/presence_histogram_group4.pdf")
plt.close()

# Abundance Count Histogram
plt.figure()
plt.xlabel("Abundance")
plt.ylabel("Samples")
plt.title("Abundance vs Samples Greater than or Equal to 600")
plt.hist(table_ged_600[:, 1])
plt.savefig(f"{figure_path}/abunance_histogram_group4.pdf")
plt.close()

# Abundance vs Presence Count Histogram
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance Greater than or Equal to 600")
plt.scatter(x=table_ged_600[:, 0], y=table_ged_600[:, 1])
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_group4.pdf")
plt.close()

"""
# Group 2
weights_leq_20 = pd.DataFrame(
    squareform(pdist(abundance_table, metric="minkowski", p=1)),
    columns=abundance_table.index,
    index=abundance_table.index)

weights = weights / np.amax(weights)
plt.figure()
plt.xlabel("Distance")
plt.ylabel("# of elements")
plt.hist(weights.to_numpy().flatten())
plt.savefig(f"{figure_path}/weight_matrix_histogram_group2.pdf")
plt.close()

# Group 3
weights = pd.DataFrame(
    squareform(pdist(abundance_table, metric="minkowski", p=1)),
    columns=abundance_table.index,
    index=abundance_table.index)

weights = weights / np.amax(weights)
plt.figure()
plt.xlabel("Distance")
plt.ylabel("# of elements")
plt.hist(weights.to_numpy().flatten())
plt.savefig(f"{figure_path}/weight_matrix_histogram_group3.pdf")
plt.close()

# Group 4
weights = pd.DataFrame(
    squareform(pdist(abundance_table, metric="minkowski", p=1)),
    columns=abundance_table.index,
    index=abundance_table.index)

weights = weights / np.amax(weights)
plt.figure()
plt.xlabel("Distance")
plt.ylabel("# of elements")
plt.hist(weights.to_numpy().flatten())
plt.savefig(f"{figure_path}/weight_matrix_histogram_group4.pdf")
plt.close()

"""
"""
# sigma = 1000000000
# kernel = np.exp(- weights ** 2 / (sigma ** 2))

kernel = weights

diagonal = pd.DataFrame(
    np.diag(np.sum(kernel)),
    columns=kernel.columns,
    index=kernel.index)

laplacian = pd.DataFrame(
    diagonal - weights,
    columns=diagonal.columns,
    index=diagonal.index)

laplacian_hat = pd.DataFrame(
    np.linalg.inv(diagonal) * laplacian,
    columns=laplacian.columns,
    index=laplacian.index)

eigenvalues, eigenvectors = np.linalg.eig(laplacian_hat)

print(eigenvalues)

plt.figure()
plt.scatter(eigenvectors[1], eigenvectors[2])
plt.savefig(f"{figure_path}/two_eigenvectors.pdf")
plt.close()
"""

# Group 2
weights = pd.DataFrame(
    squareform(pdist(abundance_table, metric="minkowski", p=1)),
    columns=abundance_table.index,
    index=abundance_table.index)

weights = weights / np.amax(weights)
plt.figure()
plt.xlabel("Distance")
plt.ylabel("# of elements")
plt.hist(weights.to_numpy().flatten(), bins=range(0, 1, 0.01))
plt.savefig(f"{figure_path}/weight_matrix_histogram_group1.pdf")
plt.close()
