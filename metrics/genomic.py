#!/usr/bin/env python3

import numpy as np
import pickle
import time
import matplotlib.pyplot as plt
import skbio

def import_otu_data(filename='../data/genetic/16S/abundance_table_97.shared', load_data=False):
    """ Import abundance data from file. sample_data is a list of strings of
    each sample name ['100001', '100002', etc.]. otu_data is a matrix of
    abundances, with otu_data[i] corresponding to sample i and consisting of a
    vector of OTUs."""

    pickle_filename = f'{filename}.pi'
    if load_data:
        with open(pickle_filename, 'rb') as f:
            sample_data, raw_otu_data = pickle.load(f)
        return sample_data, raw_otu_data

    with open(filename, 'r') as f:
        raw_otu_data = [line.strip().split("\t") for line in f][1:]
    sample_data = [val[1] for val in raw_otu_data]
    raw_otu_data = [[int(x) for x in line[3:]] for line in raw_otu_data]
    raw_otu_data = np.array(raw_otu_data)

    with open(pickle_filename, 'wb') as f:
        pickle.dump((sample_data, raw_otu_data), f)

    return sample_data, raw_otu_data

def import_metadata(filename='../data/genetic/CMAIKI_Miseq1MappingFile_20181105.csv', load_data=False):
    """ Import metadata from file. Calling metadata['Otu00001'] will return
    'Plant' or 'Fungus', etc. """

    pickle_filename = f'{filename}.pi'
    if load_data:
        with open(pickle_filename, 'rb') as f:
            metadata = pickle.load(f)
        return metadata

    with open(filename, 'r') as f:
        raw_metadata = [line.strip().split(",") for line in f][1:]
    metadata = {}
    for line in raw_metadata:
        metadata[line[0]] = line[17]

    with open(pickle_filename, 'wb') as f:
        pickle.dump(metadata, f)

    return metadata

def import_taxonomy(filename='../data/genetic/16S/annotations_97.taxonomy', load_data=False):
    """ Import OTU taxonomy from file. otu_names is a list of strings of each
    OTU ['Otu00001', 'Otu00002', 'Otu00003', etc.]. taxonomies is a dictionary
    with otu_names as the keys and values a list of taxonomic hierarchies, e.g.
    taxonomies['Otu00001'] = ['Bacteria', 'Bacterioides', etc.]. taxonomic_tree
    is a Tree class from the skbio module that contains all taxonomic
    information"""

    pickle_filename = f'{filename}.pi'
    if load_data:
        with open(pickle_filename, 'rb') as f:
            otu_names, taxonomies, taxonomic_tree  = pickle.load(f)
        return otu_names, taxonomies, taxonomic_tree

    with open(filename, 'r') as f:
        raw_metadata = [line.strip().split("\t") for line in f][1:]
    taxonomies = {}
    otu_names = []
    for line in raw_metadata:
        taxonomies[line[0]] = line[2].split(";")[:-1]
        otu_names.append(line[0])

    taxonomic_tree = skbio.TreeNode.from_taxonomy([(x, taxonomies[x]) for x in taxonomies])
    taxonomic_tree = taxonomic_tree.root_at(taxonomic_tree)
    for node in taxonomic_tree.traverse():
        node.length = 1

    #print(taxonomies['Otu00001'])
    #[print(node.name, node.distance(taxonomic_tree.root())) for node in taxonomic_tree.find('Otu00001').ancestors()]

    with open(pickle_filename, 'wb') as f:
        pickle.dump((otu_names, taxonomies, taxonomic_tree), f)
    return otu_names, taxonomies, taxonomic_tree



def get_host_distribution(sample_names, otu_data, metadata, verbose=False):
    """ Update metadata to include sample_names that weren't explicitly provided in
    the metadata file. If verbose, print the distribution of hosts. """
    hosts = {}
    count = 0
    for sample in sample_names:
        try:
            metadata[sample]
        except KeyError:
            metadata[sample] = 'no metadata'

        try:
            hosts[metadata[sample]] += 1
        except KeyError:
            hosts[metadata[sample]] = 1

    if verbose: print(hosts)
    return metadata

def plot_species_abundance_distributions(sample_names, otu_data, metadata):
    plt.figure()
    count = 0
    for i,sample in enumerate(otu_data):
        ordered_samp = np.copy(sample)
        ordered_samp.sort()
        ordered_samp = [val for val in ordered_samp[::-1] if val > 0]
        normed_samp = ordered_samp/np.sum(ordered_samp)

        #color_dict = {'Animal': 'red', 'Nonhost': 'blue', 'NA': 'orange',
        #              'Fungus': 'purple', 'Plant': 'green', 'no metadata':
        #              'None'}

        #color = color_dict[metadata[sample_names[i]]]

        if len(normed_samp) > 100:
            count += 1
            plt.step(np.linspace(0, 100, len(normed_samp)), normed_samp, lw=.5)


    lognormal = np.random.lognormal(size=500)
    ordered_samp = np.copy(lognormal)
    ordered_samp.sort()
    ordered_samp = [val for val in ordered_samp[::-1] if val > 0]
    normed_samp = ordered_samp/np.sum(ordered_samp)

    plt.step(np.linspace(0, 100, len(normed_samp)), normed_samp, color='k', lw=5)

    print(f'using {count} out of {len(otu_data)} samples')
    ax = plt.gca()
    ax.set_yscale('log')
    plt.axis([None, None, None, 1])
    plt.xlabel('\% species')
    plt.ylabel('relative abundance')
    plt.savefig('figs/species_abundance_curve_5.pdf', bbox_inches='tight')

def plot_shannon_diversities(sample_names, otu_data, metadata):
    plt.figure()
    count = 0
    shannons = []
    for i,sample in enumerate(otu_data):
        nonzero_samp = [val for val in sample if val > 0]
        normed_samp = nonzero_samp/np.sum(nonzero_samp)
        if len(normed_samp) > 100:
            shannon_div = sum([-val*np.log(val) for val in normed_samp])/np.log(len(normed_samp))
            shannons.append(shannon_div)
            count += 1

    plt.hist(shannons, bins=30, alpha=.5)

    fake_shannons = []
    for i in range(400):
        lognormal = np.random.lognormal(sigma=1.4, size=50)
        ordered_samp = [val for val in lognormal if val > 0]
        normed_samp = ordered_samp/np.sum(ordered_samp)
        shannon_div = sum([-val*np.log(val) for val in normed_samp])/np.log(len(normed_samp))
        fake_shannons.append(shannon_div)

    plt.hist(fake_shannons, bins=30, alpha=.5)

    print(f'using {count} out of {len(otu_data)} samples')
    ax = plt.gca()
    ax.set_xlim([0, 1])
    plt.xlabel('shannon diversity')
    plt.ylabel('\# occurrences')
    plt.savefig('figs/shannon_div_1.pdf', bbox_inches='tight')

def get_taxonomic_levels(taxonomic_tree):
    name_by_dist = {}
    for node in taxonomic_tree.traverse():
        try:
            name_by_dist[node.distance(taxonomic_tree.root())].append(node.name)
        except KeyError:
            name_by_dist[node.distance(taxonomic_tree.root())] = []
            name_by_dist[node.distance(taxonomic_tree.root())].append(node.name)

def get_unifrac_distances(otu_data, otu_names, taxonomic_tree, read_data=False):
    filename = 'data/unifrac_distances.pi'
    if read_data:
        with open(filename, 'rb') as f:
            unifrac_dists = pickle.load(f)
        return unifrac_dists

    unifrac_dists = skbio.diversity.beta_diversity('weighted_unifrac', otu_data,
            otu_ids=otu_names, validate=False, tree=taxonomic_tree,
            normalized=True)

    with open(filename, 'wb') as f:
        pickle.dump(unifrac_dists, f)

    return unifrac_dists

def get_unweighted_unifrac_distances(otu_data, otu_names, taxonomic_tree, read_data=False):
    filename = 'data/unweighted_unifrac_distances.pi'
    if read_data:
        with open(filename, 'rb') as f:
            unifrac_dists = pickle.load(f)
        return unifrac_dists

    unifrac_dists = skbio.diversity.beta_diversity('unweighted_unifrac', otu_data,
            otu_ids=otu_names, validate=False, tree=taxonomic_tree)

    with open(filename, 'wb') as f:
        pickle.dump(unifrac_dists, f)

    return unifrac_dists

def get_l2_distances(otu_data, otu_names, taxonomic_tree, read_data=False):
    filename = 'data/l2_distances.pi'
    if read_data:
        with open(filename, 'rb') as f:
            l2_dists = pickle.load(f)
        return l2_dists

    l2_dists = []
    for i,samp_1 in enumerate(otu_data):
        print(i)
        for j,samp_2 in enumerate(otu_data):
            l2_dists.append(np.log(np.linalg.norm(samp_1-samp_2)))

    l2_dists = np.array(l2_dists)
    l2_dists = np.reshape(l2_dists, (len(otu_data), len(otu_data)))

    with open(filename, 'wb') as f:
        pickle.dump(l2_dists, f)

    return l2_dists

def get_dists(otu_data, otu_names, taxonomic_tree):
    unifrac_dists = get_unifrac_distances(otu_data, otu_names, taxonomic_tree,
            read_data=True)

    unweighted_unifrac_dists = get_unweighted_unifrac_distances(otu_data, otu_names, taxonomic_tree,
            read_data=True)

    l2_dists = get_l2_distances(otu_data, otu_names, taxonomic_tree,
            read_data=False)

    return unifrac_dists.redundant_form(), unweighted_unifrac_dists.redundant_form(), l2_dists

def plot_heatmaps(unifrac_dists, unweighted_unifrac_dists, l2_dists):
    for name, dist in zip(['unifrac', 'unweighted', 'l2'], [unifrac_dists,
                            unweighted_unifrac_dists, l2_dists]):
        plt.figure()
        plt.imshow(dist, cmap='hot')
        plt.xlabel('sample \#')
        plt.ylabel('sample \#')
        cbar = plt.colorbar()
        cbar.set_label('distance')
        plt.title(f'{name}')
        plt.savefig(f'figs/{name}_distance_matrix.pdf', bbox_inches='tight')

if __name__ == '__main__':
    sample_names, otu_data = import_otu_data(load_data=False)
    metadata = import_metadata(load_data=False)
    metadata = get_host_distribution(sample_names, otu_data, metadata, verbose=True)
    otu_names, taxonomies, taxonomic_tree = import_taxonomy(load_data=False)

    plot_species_abundance_distributions(sample_names, otu_data, metadata)
    plot_shannon_diversities(sample_names, otu_data, metadata)

    unifrac_dists, unweighted_unifrac_dists, l2_dists = get_dists(otu_data, otu_names, taxonomic_tree)
    plot_heatmaps(unifrac_dists, unweighted_unifrac_dists, l2_dists)
    node_by_dist = get_taxonomic_levels(taxonomic_tree)
    
