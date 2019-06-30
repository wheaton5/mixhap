#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description = "haplotype phasing via sparse mixture model clustering")

parser.add_argument("-i","--variants", required=True, help = "variant read support file")
parser.add_argument("-w","--window", required=False, type=int, default=100, help = "window size of number of variants to consider")
parser.add_argument("-s","--stepsize", required=False, type=int, default=10, help = "step size to move the window each iteration")
parser.add_argument("-r","--retries", required=False, type=int, default=10, help = "number of retries of cluster center initialization (only affects first window)")


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
            
    variant_data_new = []
    for i in range(step_size):
        index = i + len(variant_data_so_far)
        try:
            toks = input_data.readline().strip().split()
            molecule_data_new = []
            for tok in toks:
                
                subtok = tok.split(":")
                molecule_id = int(tok[0])
                allele_fraction = float(tok[1])
                if not molecule_id in molecule_indexes:
                    molecule_indexes[molecule_id] = len(molecule_variant_indices)
                    molecule_variant_indices.append([])
                    molecule_variant_allele_fractions.append([])
                molecule_index = molecule_indexes[molecule_id]
                molecule_variant_indices[molecule_index].append(index)
                molecule_variant_allele_fractions[molecule_index].append(allele_fraction)
                molecule_data_new.append((molecule_id, allele_fraction))
            variant_data_new.append(molecule_data_new)
    return (molecule_variant_indices, molecule_variant_allele_fractions, variant_data_so_far) 
                



input_data = open(args.variants)
variant_data_so_far = [] # should be a list of variants with each variant being a list of (molecule_id, allele_fraction)
molecule_variant_indices = []
molecule_variant_allele_fractions = np.zeros(0)
molecule_ids = []
while True:
    molecule_variant_indices, molecule_variant_fractions, variant_data_so_far = 
            collect_data(input_data, variant_data_so_far, args.window, args.step_size)
    if len(molecules) == 0:
        break
