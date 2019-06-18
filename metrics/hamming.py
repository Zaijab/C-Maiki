#!/usr/bin/env python
import numpy as np

f = open("../data/genetic/16S/sequences_97.fasta")
all_otu_data = f.read()

otu_list = all_otu_data.split(">")

otu_dna_pair = []

for otu in otu_list:
    otu_dna_pair.append(otu.split("\n"))

def hamming_dist(otu1, otu2):
    count = 0
    for i in range(len(otu1)):
        if otu1[i] != otu2[i]:
            count += 1
    return count


def calc_dist_matrix(otu_dna_pair_list):
    dist_matrix = [[0 for x in range(len(otu_dna_pair_list))] for y in range(len(otu_dna_pair_list))]
    for i in range(len(otu_dna_pair_list)):
        for j in range(len(otu_dna_pair_list)):
            dist_matrix[i][j] = hamming_dist(otu_dna_pair_list[i][1],otu_dna_pair_list[j][1])
            #print(dist_matrix[i][j])
    return dist_matrix

hamming_dist_matrix = calc_dist_matrix(otu_dna_pair[1:])

np.save("hamming_dist_matrix",hamming_dist_matrix)
print(hamming_dist_matrix)
