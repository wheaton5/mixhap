#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description = "haplotype phasing via sparse mixture model clustering")

parser.add_argument("-i","--variants", required=True, help = "variant read support file")
parser.add_argument("-w","--window", required=False, type=int, default = 100, help = "window size of number of variants to consider")
parser.add_argument("-s","--step_size", required=False, type=int, default = 10, help = "step size to move the window each iteration")
parser.add_argument("-r","--retries", required=False, type=int, default = 10, help = "number of retries of cluster center initialization (only affects first window)")
parser.add_argument("-p","--ploidy", required=False, type=int, default = 2, help = "ploidy") 
parser.add_argument("--iterations", required=False, type=int, default=10, help="iterations of sparse fuzzy kmeans algorithm")
args = parser.parse_args()

import numpy as np
from scipy.special import logsumexp

# build the data we need for next round
def collect_data(input_data, variant_data_so_far, window, step_size):
    molecule_indexes = {}
    molecule_variant_indices = []
    molecule_variant_allele_fractions = []
     
    if len(variant_data_so_far) >= window:
        variant_data_so_far = variant_data_so_far[step_size-1:]
    for (index, variant) in enumerate(variant_data_so_far):
        for (molecule_id, allele_fraction) in variant:
            if not molecule_id in molecule_indexes:
                molecule_indexes[molecule_id] = len(molecule_variant_indices)
                molecule_variant_indices.append([])
                molecule_variant_allele_fractions.append([])
            molecule_index = molecule_indexes[molecule_id]
            molecule_variant_indices[molecule_index].append(index)
            molecule_variant_allele_fractions[molecule_index].append(allele_fraction)
    var_offset = len(variant_data_so_far)
    for i in range(step_size):
        index = i + var_offset
        line = input_data.readline()
        if line == '':
            if i > 0:
                return (molecule_variant_indices, molecule_variant_allele_fractions, variant_data_so_far, False)
            else:
                return (None, None, None, True)
        toks = line.strip().split()
        molecule_data_new = []
        for tok in toks:
            subtok = tok.split(":")
            molecule_id = int(subtok[0])
            allele_fraction = float(subtok[1])
            if not molecule_id in molecule_indexes:
                molecule_indexes[molecule_id] = len(molecule_variant_indices)
                molecule_variant_indices.append([])
                molecule_variant_allele_fractions.append([])
            molecule_index = molecule_indexes[molecule_id]
            molecule_variant_indices[molecule_index].append(index)
            molecule_variant_allele_fractions[molecule_index].append(allele_fraction)
            molecule_data_new.append((molecule_id, allele_fraction))
        variant_data_so_far.append(molecule_data_new)
    for index in range(len(molecule_variant_allele_fractions)):
        molecule_variant_allele_fractions = np.matrix(molecule_variant_allele_fractions)
    return (molecule_variant_indices, molecule_variant_allele_fractions, variant_data_so_far, False) 

input_data = open(args.variants)
variant_data_so_far = [] # should be a list of variants with each variant being a list of (molecule_id, allele_fraction)
np.random.seed(4) # guaranteed random seed chosen by dice roll (joke) https://xkcd.com/221/
haplotype_centers_sum = np.zeros((args.ploidy, args.step_size))
haplotype_centers_denom = np.zeros((args.ploidy, args.step_size)) + 1e-7 # epsilon to avoid div by 0
data_grabs = 0
previous = False
while True: # probably eventually loop over connected components of the graph, or over chromosomes
    (molecule_variant_indices, molecule_variant_fractions, variant_data_so_far, stop) = \
            collect_data(input_data, variant_data_so_far, args.window, args.step_size)
    data_grabs += 1
    print("data grabs",data_grabs)
    if stop:
        break
    
    best_centers = None#haplotype_centers.copy()
    best_centers_distance = 10000000000
    for repitition in range(args.retries):
        if previous:
            haplotype_centers = np.concatenate((haplotype_centers[:,min(haplotype_centers.shape[1]-1), args.step_size-1:].T, 
                np.tile([0.5], (args.step_size, args.ploidy)) + np.random.rand(args.ploidy, args.step_size) / 20.0 )) 
        else:
            haplotype_centers = np.random.rand(args.ploidy, args.step_size)
            previous = True
        #print("init",haplotype_centers)
        for iteration in range(args.iterations):
            center_distance = 0
            for (mol_indices, mol_fracs) in zip(molecule_variant_indices, molecule_variant_fractions):
                mol_distances = np.sum(np.square(mol_fracs - haplotype_centers[:,mol_indices]), axis=1) # standard sum of square differences
                mol_distances /= np.sum(mol_distances) # normalize to sum to 1, now this can be used as a probability for this molecule to each cluster
                mol_distances = 1.0 - mol_distances
                mol_distances = np.matrix(mol_distances) # more useful shape for this object
                haplotype_centers_sum[:,mol_indices] += mol_distances * mol_fracs # matrix multiplication to produce weighted sum
                haplotype_centers_denom[:,mol_indices] += np.broadcast_to(mol_distances, (args.ploidy, len(mol_indices))) # and the weights for later average
                center_distance += logsumexp(np.log(mol_distances))
            haplotype_centers = np.divide(haplotype_centers_sum , haplotype_centers_denom)
            haplotype_centers_sum = np.zeros((args.ploidy, args.step_size))
            haplotype_centers_denom = np.zeros((args.ploidy, args.step_size)) + 1e-7 # epsilon to avoid div by 0
            
        if center_distance < best_centers_distance:
            best_centers = haplotype_centers.copy()
            best_centers_distance = center_distance
    print("best centers",best_centers)


            
            
