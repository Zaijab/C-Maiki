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

    # Plotting Functions

    def plot_heatmaps(self):
        for name, dist in zip(['unifrac', 'unweighted', 'l2'], [self.unifrac_dists, self.unweighted_unifrac_dists, self.l2_dists]):
            plt.figure()
            plt.imshow(dist, cmap='hot')
            plt.xlabel('sample \#')
            plt.ylabel('sample \#')
            cbar = plt.colorbar()
            cbar.set_label('distance')
            plt.title(f'{name}')
            plt.savefig(
                f'{self.figure_path}/{name}_distance_matrix.pdf', bbox_inches='tight')

    def plot_shannon_diversities(self):
        plt.figure()
        count = 0
        shannons = []
        for i, sample in enumerate(self.otu_data):
            nonzero_samp = [val for val in sample if val > 0]
            normed_samp = nonzero_samp/np.sum(nonzero_samp)
            if len(normed_samp) > 100:
                shannon_div = sum([-val*np.log(val)
                                   for val in normed_samp])/np.log(len(normed_samp))
                shannons.append(shannon_div)
                count += 1
        plt.hist(shannons, bins=30, alpha=.5)
        fake_shannons = []
        for i in range(400):
            lognormal = np.random.lognormal(sigma=1.4, size=50)
            ordered_samp = [val for val in lognormal if val > 0]
            normed_samp = ordered_samp/np.sum(ordered_samp)
            shannon_div = sum([-val*np.log(val)
                               for val in normed_samp])/np.log(len(normed_samp))
            fake_shannons.append(shannon_div)
        plt.hist(fake_shannons, bins=30, alpha=.5)
        print(f'using {count} out of {len(self.otu_data)} samples')
        ax = plt.gca()
        ax.set_xlim([0, 1])
        plt.xlabel('shannon diversity')
        plt.ylabel('\# occurrences')
        plt.savefig(f'{self.figure_path}/shannon_div_1.pdf',
                    bbox_inches='tight')

    def plot_species_abundance_distributions(self):
        plt.figure()
        count = 0
        for i, sample in enumerate(self.otu_data):
            ordered_samp = np.copy(sample)
            ordered_samp.sort()
            ordered_samp = [val for val in ordered_samp[::-1] if val > 0]
            normed_samp = ordered_samp/np.sum(ordered_samp)
            if len(normed_samp) > 100:
                count += 1
                plt.step(np.linspace(0, 100, len(normed_samp)),
                         normed_samp, lw=.5)
        lognormal = np.random.lognormal(size=500)
        ordered_samp = np.copy(lognormal)
        ordered_samp.sort()
        ordered_samp = [val for val in ordered_samp[::-1] if val > 0]
        normed_samp = ordered_samp/np.sum(ordered_samp)
        plt.step(np.linspace(0, 100, len(normed_samp)),
                 normed_samp, color='k', lw=5)
        print(f'using {count} out of {len(self.otu_data)} samples')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.axis([None, None, None, 1])
        plt.xlabel('\% species')
        plt.ylabel('relative abundance')
        plt.savefig(
            f'{self.data_path}/species_abundance_curve_5.pdf', bbox_inches='tight')

    def plot_abundance_histogram(self):
        plt.figure()
        # Otu Count - otu_count[i] = total number of otu in all sample
        self.abundance_count = np.sum(self.otu_data, axis=1)
        #axes = plt.gca()
        #axes.set_xlim([10, 10000])
        plt.xlabel("Abundance")
        plt.ylabel("Sample")
        plt.hist(self.abundance_count, bins=100)
        plt.savefig(f'{self.figure_path}/abundance_histogram_sana.pdf')
        plt.close()

    def plot_presence_histogram(self):
        plt.figure()
        # Otu Count - otu_count[i] = total number of otu in all sample
        self.presence_count = np.sum(
            np.where(self.otu_data > 0, 1, self.otu_data), axis=1)
        axes = plt.gca()
        axes.set_xlim([20, 10000])
        plt.xlabel("Presence")
        plt.ylabel("Sample")
        plt.hist(self.presence_count, bins=100)
        plt.savefig(f'{self.figure_path}/presence_histogram_sana.pdf')
        plt.close()

    def plot_presence_vs_abundance_scatter(self):
        plt.figure()
        plt.xlabel("Presence")
        plt.ylabel("Abundance")
        plt.scatter(self.presence_count, self.abundance_count)
        plt.savefig(
            f'{self.figure_path}/presence_vs_abundance_scatter_sana.pdf')
        plt.close()


def main():
    data_path = "./amandin/data/microbiome/sana"
    figure_path = "./amandin/presentation/microbiome/figures"
    cache_path = "./amandin/data/microbiome/cache"
    load_data = False
    test_analysis = Microbiome(data_path, figure_path, cache_path, load_data)

    test_analysis.plot_abundance_histogram()
    test_analysis.plot_presence_histogram()
    test_analysis.plot_presence_vs_abundance_scatter()


main()
