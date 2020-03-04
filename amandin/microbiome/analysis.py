# Import dependencies
import h5py
import numpy as np
import pickle
import time
import matplotlib.pyplot as plt

# Initialize Files and Paths
microbiome = h5py.File("micorbiome.hdf5", "w")
otu_data_path = "../data/microbiome/zain/abundance_table_97.shared"

# Import Relevent Data from files


def import_otu_data(filename='../data/microbiome/zain/abundance_table_97.shared', load_data=False):
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
            otu_names, taxonomies, taxonomic_tree = pickle.load(f)
        return otu_names, taxonomies, taxonomic_tree

    with open(filename, 'r') as f:
        raw_metadata = [line.strip().split("\t") for line in f][1:]
    taxonomies = {}
    otu_names = []
    for line in raw_metadata:
        taxonomies[line[0]] = line[2].split(";")[:-1]
        otu_names.append(line[0])

    taxonomic_tree = skbio.TreeNode.from_taxonomy(
        [(x, taxonomies[x]) for x in taxonomies])
    taxonomic_tree = taxonomic_tree.root_at(taxonomic_tree)
    for node in taxonomic_tree.traverse():
        node.length = 1
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

    if verbose:
        print(hosts)
    return metadata


# Import Otu Data and metadata
data_path = "./amandin/data/microbiome/zain/abundance_table_100.shared"
sample_names, otu_data = import_otu_data(data_path)

# Start creating Abundance table with metadata
microbiome = h5py.File("mytestfile.hdf5", "w")
microbiome.create_dataset("otu_data", otu_data)
microbiome.close()
