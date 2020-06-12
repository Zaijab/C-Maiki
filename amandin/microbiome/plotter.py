
# Filtering
diverse = abundance_table.loc[500 <= abundance_table["Presence"]]
diverse_metadata = metadata.loc[diverse.index]

diverse_metadata.to_csv(f"{data_path}/presence_above_500.csv")

site_name = list(host_plant_group_metadata["site_name"])

for i, site in enumerate(site_name):
    if "Summit" in site:
        site_name[i] = 1
    if "Gardens" in site:
        site_name[i] = 2
    if "STFalls" in site:
        site_name[i] = 3
    if "Estuary" in site:
        site_name[i] = 4
    if "DrumRoad" in site:
        site_name[i] = 5
    if "AboveFalls" in site:
        site_name[i] = 6

plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance Plantroot")
plt.scatter(host_plant_abundance_table["Presence"],
            host_plant_abundance_table["Abundance"], label='Plantroot', c=site_name)
plt.savefig(
    f"{figure_path}/presence_vs_abundance_scatter_plantroot_site_name_color.png")
plt.close()

plt.figure()
plt.xlabel("Presence")
plt.ylabel("Sample")
plt.title("Presence Histogram")
plt.hist(plantroot_group["Presence"], bins=80)
plt.savefig(f"{figure_path}/plantroot_presence_histogram_80_bins.png")
plt.close()

normalized_group_weight = [
    weight / np.amax(weight.to_numpy()) for weight in group_weight]
"""
groups = [
    abundance_table,
    abundance_table.loc[abundance_table["Presence"] < 20],
    abundance_table.loc[(20 <= abundance_table["Presence"])
                        & (abundance_table["Presence"] < 600)],
    abundance_table.loc[600 <= abundance_table["Presence"]],
]
"""

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
abundance_table.drop(
    'Abundance', axis=1

for sigma in np.linspace(0, 1, 100):
    for i, dists in enumerate(group_weights):
        dist=dists / np.amax(dists)
        dist=np.exp(- (dist ** 2) / (sigma ** 2))
        plt.figure()
        plt.xlabel("Kernel")
        plt.ylabel("# of Elements")
        plt.title(f"Kernel Group {i+1} Sigma {sigma}")
        plt.hist(dist.to_numpy().flatten(), bins=np.linspace(0, 1, 50))
        plt.savefig(
            f"{figure_path}/sigmas/kernel_matrix_histogram_group{i+1}_{sigma}.pdf")
        plt.close()
"""
for i, distance_matrix in enumerate(group_weights):
    dist = distance_matrix / np.amax(distance_matrix)
    mu = np.mean(dist)
    variance = np.var(dist)
    dist = ((dist - mu) ** 2) / variance
    plt.figure()
    plt.xlabel("Kernel")
    plt.ylabel("# of Elements")
    plt.title(f"Kernel Group {i+1} Sigma {sigma}")
    plt.hist(dist.to_numpy().flatten(), bins=np.linspace(0, 1, 50))
    plt.savefig(
        f"{figure_path}/sigmas/kernel_matrix_histogram_group{i+1}_zscore.pdf")
    plt.close()

"""
"""
# Only look at Plantroot groups
plantroot_ids = metadata.loc[metadata["collection_label"] == "PlantRoot"]
plantroot_group = abundance_table.filter(
    items=list(plantroot_ids.index), axis=0)
plantroot_group_metadata = plantroot_ids.loc[plantroot_group.index]
groups = [plantroot_group]
plantroot_group_metadata.to_csv(f"{data_path}/plantroot_group_metadata.csv")


# Computation Zone
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
    weight / 14835 for weight in group_weight]

for weight in normalized_group_weight:
    plt.figure()
    plt.xlabel("Distance")
    plt.ylabel("Sample")
    plt.title("Distance Histogram")
    A = weight.to_numpy()
    plt.hist(A[~np.eye(A.shape[0], dtype=bool)].reshape(
        A.shape[0], -1).flatten(), bins=80)
    plt.savefig(f"{figure_path}/plantroot_normalized_distance_histogram.png")
    plt.close()

# Compute Kernel exp(-(w_ij)^2/sigma^2)
group_kernel = [np.exp(- (weight - np.mean(weight.to_numpy())) ** 2 / np.var(weight.to_numpy()))
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

# Compute Normalized Laplacian
normalized_group_laplacian = [
    pd.DataFrame(
        np.linalg.inv(diagonal) * laplacian,
        columns=diagonal.columns,
        index=diagonal.index)
    for diagonal, laplacian in zip(group_diagonal, group_laplacian)]

# Compute Eigenvectors / Eigenvalues
group_eigen = [eigs(laplacian.to_numpy(), k=4)
               for laplacian in group_laplacian]

print(group_eigen[0][1].imag.T)
# Plot
for i, eigen in enumerate(group_eigen):
    eigenvalues, eigenvectors = eigen[0].real, eigen[1].real.T
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(eigenvectors[1], eigenvectors[2], eigenvectors[3])
    plt.title("Laplacian Eigs Plantroot 3d")
    for label, x, y, z in zip(plantroot_group_metadata["site_code"], eigenvectors[1], eigenvectors[2], eigenvectors[3]):
        ax.text(x, y, z, label, None)
    plt.savefig(
        f"{figure_path}/two_eigenvectors_non_normalized_eigs_plantroot3d.png")
    plt.close()
"""
plt.figure()
plt.scatter(eigenvectors[1], eigenvectors[2])
plt.title("Laplacian Eigs Mosquito")
plt.savefig(
    f"{figure_path}/two_eigenvectors_non_normalized_eigs_mosquito_2d.png")
plt.close()
# Host - Plant
host_plant_ids=metadata.loc[metadata["host"] == "Plant"]
host_plant_abundance_table=abundance_table.filter(
    items=list(host_plant_ids.index), axis=0)
host_plant_group_metadata=host_plant_ids.loc[host_plant_abundance_table.index]

# Plotting
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples Host Plant")
plt.hist(plant_group["Presence"], bins=100)
plt.savefig(f"{figure_path}/presence_histogram_host_plant.png")
plt.close()

plt.figure()
plt.xlabel("Abundance")
plt.ylabel("Samples")
plt.title("Abundance vs Samples Host Nonhost")
plt.hist(host_nonhost_abundance_table["Abundance"], bins=100)
plt.savefig(f"{figure_path}/abundance_histogram_host_nonhost.png")
plt.close()

plt.figure()
plt.xlabel("Presence")
plt.ylabel("Samples")
plt.title("Presence vs Samples Host Nonhost")
plt.hist(host_nonhost_abundance_table["Presence"], bins=100)
plt.savefig(f"{figure_path}/presence_histogram_host_nonhost.png")
plt.close()

plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance Host Nonhost")
plt.scatter(host_nonhost_abundance_table["Presence"],
            host_nonhost_abundance_table["Abundance"])
plt.savefig(
    f"{figure_path}/presence_vs_abundance_scatter_host_nonhost.png")
plt.close()


# Plotting
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance")
plt.scatter(host_nonhost_abundance_table["Presence"],
            host_nonhost_abundance_table["Abundance"], alpha=0.5, label='Nonhost')
plt.scatter(host_plant_abundance_table["Presence"],
            host_plant_abundance_table["Abundance"], alpha=0.5, label='Plant', color='green')
plt.scatter(host_animal_abundance_table["Presence"],
            host_animal_abundance_table["Abundance"], alpha=0.5, label='Animal')
plt.legend(loc='upper right')
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_hosts.png")
plt.close()


# Hosts
host_nonhost_ids=metadata.loc[metadata["host"] == "Nonhost"]
host_nonhost_abundance_table=abundance_table.filter(
    items=list(host_nonhost_ids.index), axis=0)
print("Nonhost:", host_nonhost_abundance_table.shape)

host_plant_ids=metadata.loc[metadata["host"] == "Plant"]
host_plant_abundance_table=abundance_table.filter(
    items=list(host_plant_ids.index), axis=0)
print("Plant:", host_plant_abundance_table.shape)

host_animal_ids=metadata.loc[metadata["host"] == "Animal"]
host_animal_abundance_table=abundance_table.filter(
    items=list(host_animal_ids.index), axis=0)
print("Animal:", host_animal_abundance_table.shape)

# Plotting
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance")
plt.scatter(host_nonhost_abundance_table["Presence"],
            host_nonhost_abundance_table["Abundance"], alpha=0.5, label='Nonhost')
plt.scatter(host_plant_abundance_table["Presence"],
            host_plant_abundance_table["Abundance"], alpha=0.5, label='Plant', color='green')
plt.scatter(host_animal_abundance_table["Presence"],
            host_animal_abundance_table["Abundance"], alpha=0.5, label='Animal')
plt.legend(loc='upper right')
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_hosts.png")
plt.close()

host_plant_ids=metadata.loc[metadata["collection_label"] == "PlantRoot"]
host_plant_abundance_table=abundance_table.filter(
    items=list(host_plant_ids.index), axis=0)

host_nonplant_ids=metadata.loc[metadata["collection_label"] != "PlantRoot"]
host_nonplant_abundance_table=abundance_table.filter(
    items=list(host_nonplant_ids.index), axis=0)

plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance")
plt.scatter(host_nonplant_abundance_table["Presence"],
            host_nonplant_abundance_table["Abundance"], alpha=1, label='Not Plantroot')
plt.scatter(host_plant_abundance_table["Presence"],
            host_plant_abundance_table["Abundance"], alpha=1, label='Plant', color='green')
plt.legend(loc='upper right')
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_plantroot.png")
plt.close()

# Habitats
habitat_marine_ids=metadata.loc[metadata["habitat"] == "Marine"]
habitat_marine_abundance_table=abundance_table.filter(
    items=list(habitat_marine_ids.index), axis=0)
print("Marine:", habitat_marine_abundance_table.shape)

habitat_terrestrial_ids=metadata.loc[metadata["habitat"] == "Terrestrial"]
habitat_terrestrial_abundance_table=abundance_table.filter(
    items=list(habitat_terrestrial_ids.index), axis=0)
print("Terrestrial:", habitat_terrestrial_abundance_table.shape)

habitat_riverine_ids=metadata.loc[metadata["habitat"] == "Riverine"]
habitat_riverine_abundance_table=abundance_table.filter(
    items=list(habitat_riverine_ids.index), axis=0)
print("Riverine:", habitat_riverine_abundance_table.shape)

# Plotting
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Sample")
plt.title("Presence vs Sample")
plt.hist(
    habitat_Marine_abundance_table["Presence"], bins=100, alpha=0.5, label='Marine')
plt.hist(habitat_terrestrial_abundance_table["Presence"],
         bins=100, alpha=0.5, label='Terrestrial', color='green')
plt.hist(habitat_riverine_abundance_table["Presence"],
         bins=100, alpha=0.5, label='Riverine')
plt.legend(loc='upper right')
plt.savefig(f"{figure_path}/presence_histogram_habitats.png")
plt.close()

# Plotting
plt.figure()
plt.xlabel("Presence")
plt.ylabel("Abundance")
plt.title("Presence vs Abundance")
plt.hist(habitat_Marine_abundance_table["Presence"],
         habitat_Marine_abundance_table["Abundance"], alpha=0.5, label='Marine')
plt.hist(habitat_terrestrial_abundance_table["Presence"],
         habitat_terrestrial_abundance_table["Abundance"], alpha=0.5, label='Terrestrial', color='green')
plt.hist(habitat_riverine_abundance_table["Presence"],
         habitat_riverine_abundance_table["Abundance"], alpha=0.5, label='Riverine')
plt.legend(loc='upper right')
plt.savefig(f"{figure_path}/presence_vs_abundance_scatter_habitats.png")
plt.close()

diagonal=host_plant_group_metadata["lat"] + host_plant_group_metadata["long"]
low_right=min(diagonal)
high_left=max(diagonal)
gps_delta=high_left - low_right
diagonal=(diagonal - low_right) / gps_delta
color_gradient=[(0, i, 0) for i in list(diagonal)]

stfalls_mosqu=mosquito_metadata.loc[mosquito_metadata["transect_name"] == "STFalls"]
stfalls_abund=abundance_table.filter(
    items=list(stfalls_mosqu.index), axis=0)
