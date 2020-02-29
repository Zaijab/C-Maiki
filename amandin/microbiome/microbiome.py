"""
This module is created to do data analysis on microbiome data.
"""

from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import pickle
import time
import skbio


class Microbiome(object):
    """
    Represents a microbiome
    Used to input microbial data and perform analysis on it.
    """

    # Initializing Functions

    def __init__(self, data_path, figure_path, cache_path, load_data):
        """
        Construct a new microme
        """
        self.data_path = data_path
        self.figure_path = figure_path
        self.cache_path = cache_path
        self.load_data = load_data
        self.sample_data, self.otu_data = self.import_otu_data()
        self.metadata = self.import_metadata()
        self.otu_names, self.taxonomies, self.taxonomic_tree = self.import_taxonomy()

    def import_otu_data(self):
        """
        Import Abundance Data from data_path.

        Sample_data - a list of strings of
        each sample name ['100001', '100002', etc.].

        otu_data - a matrix of abundances
        otu_data[i] corresponding to sample i and consisting of a
        vector of OTUs.
        """

        # Use filename
        otu_path = f"{self.data_path}/abundance_table_97.shared"

        # Import procedure
        with open(otu_path, 'r') as f:
            raw_otu_data = [line.strip().split("\t") for line in f][1:]
        sample_data = [val[1] for val in raw_otu_data]
        raw_otu_data = [[int(x) for x in line[3:]] for line in raw_otu_data]
        raw_otu_data = np.array(raw_otu_data)

        return sample_data, raw_otu_data

    def import_metadata(self):
        """
        Import metadata from file.
        Calling metadata['Otu00001'] will return 'Plant' or 'Fungus', etc.
        """

        # Use filename
        metadata_path = f'{self.data_path}/SuperTransect_mapping_file.csv'

        # Import Procedure
        with open(metadata_path, 'r') as f:
            raw_metadata = [line.strip().split(",") for line in f][1:]
        metadata = {}
        for line in raw_metadata:
            metadata[line[0]] = line[17]

        hosts = {}
        count = 0
        for sample in self.sample_data:
            try:
                metadata[sample]
            except KeyError:
                metadata[sample] = 'no metadata'

            try:
                hosts[metadata[sample]] += 1
            except KeyError:
                hosts[metadata[sample]] = 1
        return metadata

    def import_taxonomy(self):
        """
        Import OTU taxonomy from file.

        Output:
        otu_names = ['Otu00001', 'Otu00002', 'Otu00003', etc.].
        taxonomies['Otu00001'] = ['Bacteria', 'Bacterioides', etc.].
        taxonomic_tree is a Tree class from the skbio module that contains all taxonomic information
        """

        # Use filename
        taxonomy_path = f"{self.data_path}/annotations_97.taxonomy"

        # Import Procedure
        with open(taxonomy_path, 'r') as f:
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

        return otu_names, taxonomies, taxonomic_tree

    # Analysis Function

    def get_unifrac_distances(self):
        unifrac_dists = skbio.diversity.beta_diversity(
            'weighted_unifrac', self.otu_data, otu_ids=self.otu_names, validate=False, tree=self.taxonomic_tree, normalized=True)
        return unifrac_dists

    def get_unweighted_unifrac_distances(self):
        unifrac_dists = skbio.diversity.beta_diversity(
            'unweighted_unifrac', self.otu_data, otu_ids=self.otu_names, validate=False, tree=self.taxonomic_tree)

        return unifrac_dists

    def get_l2_distances(self):
        l2_dists = []
        for i, samp_1 in enumerate(self.otu_data):
            print(i)
            for j, samp_2 in enumerate(self.otu_data):
                l2_dists.append(np.log(np.linalg.norm(samp_1-samp_2)))
        l2_dists = np.array(l2_dists)
        l2_dists = np.reshape(
            l2_dists, (len(self.otu_data), len(self.otu_data)))

        return l2_dists

    def get_dists(self):
        unifrac_dists = self.get_unifrac_distances()
        unweighted_unifrac_dists = self.get_unweighted_unifrac_distances()
        l2_dists = self.get_l2_distances()
        self.unifrac_dists = unifrac_dists.redundant_form()
        self.unweighted_unifrac_dists = unweighted_unifrac_dists.redundant_form()
        self.l2_dists = l2_dists
        return unifrac_dists.redundant_form(), unweighted_unifrac_dists.redundant_form(), l2_dists

    def l1_distances(self, vals, name):
        plt.figure()
        plt.title("Heatmap l1 distance")
        plt.ylabel("Sample")
        plt.xlabel("Distance")
        dists = distance.pdist(vals, metric='minkowski', p=1)
        self.dists = dists
        plt.hist(dists)
        plt.savefig(f'{self.figure_path}/l1_dist_{name}.pdf')

    def grouper(self):
        # Create 3 groups based on presence
        left, right = 20, 800
        groups = [[] for i in range(3)]
        for i in self.otu_data:
            precence = np.sum(np.where(i > 0, 1, i))
            if precence < left:
                groups[0].append(i)
            elif left < precence < right:
                groups[1].append(i)
            else:
                groups[2].append(i)
        return [groups[0], groups[1], groups[2]]

    def compute_abundance(self, vals):
        self.abunance = np.sum(vals, axis=1)
        return self.abunance

    def compute_precence(self, vals):
        self.precence = np.sum(np.where(vals > 0, 1, vals), axis=1)
        return self.precence
    # Plotting Functions

    def plot_abundance_histogram(self, vals, name=None):
        plt.figure()
        #axes = plt.gca()
        #axes.set_xlim([10, 10000])
        plt.xlabel("Abundance")
        plt.ylabel("Sample")
        plt.hist(np.sum(vals, axis=1), bins=100)
        plt.savefig(f'{self.figure_path}/abundance_histogram_{name}.pdf')
        plt.close()

    def plot_presence_histogram(self, vals, name=None):
        plt.figure()
        #axes = plt.gca()
        #axes.set_xlim([20, 10000])
        plt.xlabel("Presence")
        plt.ylabel("Sample")
        plt.hist(np.sum(np.where(vals > 0, 1, vals), axis=1), bins=100)
        plt.savefig(f'{self.figure_path}/presence_histogram_{name}.pdf')
        plt.close()

    def plot_presence_vs_abundance_scatter(self, vals, name=None):
        plt.figure()
        plt.xlabel("Presence")
        plt.ylabel("Abundance")
        #axes = plt.gca()
        #axes.set_xlim([-1, 25])
        abundance = np.sum(vals, axis=1)
        precence = np.sum(np.where(vals > 0, 1, vals), axis=1)
        plt.scatter(precence, abundance)
        plt.savefig(
            f'{self.figure_path}/presence_vs_abundance_scatter_{name}.pdf')
        plt.close()


def main():
    data_path = "./amandin/data/microbiome/zain"
    figure_path = "./amandin/presentation/microbiome/figures"
    cache_path = "./amandin/data/microbiome/cache"
    load_data = False
    test_analysis = Microbiome(data_path, figure_path, cache_path, load_data)

    otu_data = test_analysis.otu_data


"""
    test_analysis.plot_abundance_histogram(otu_data, "total")
    test_analysis.plot_presence_histogram(otu_data, "total")
    test_analysis.plot_presence_vs_abundance_scatter(otu_data, "total")

    groups = test_analysis.grouper()

    test_analysis.plot_presence_histogram(np.array(groups[0]), "less_than_20")
    test_analysis.plot_abundance_histogram(np.array(groups[0]), "less_than_20")
    test_analysis.plot_presence_vs_abundance_scatter(
        np.array(groups[0]), "less_than_20")

    test_analysis.plot_presence_histogram(np.array(groups[1]), "between")
    test_analysis.plot_abundance_histogram(np.array(groups[1]), "between")
    test_analysis.plot_presence_vs_abundance_scatter(
        np.array(groups[1]), "between")

    print(type(groups[2]))
    test_analysis.plot_presence_histogram(np.array(groups[2]), "more_than_800")
    test_analysis.plot_abundance_histogram(
        np.array(groups[2]), "more_than_800")
    test_analysis.plot_presence_vs_abundance_scatter(
        np.array(groups[2]), "more_than_800")
"""

test_analysis.l1_distances(groups[0], "Less than 20")
test_analysis.l1_distances(groups[1], "Between")
test_analysis.l1_distances(groups[2], "Greater than 800")

print(len(groups[1]), len(groups[1])*(len(groups[1])-1)/2)


main()
