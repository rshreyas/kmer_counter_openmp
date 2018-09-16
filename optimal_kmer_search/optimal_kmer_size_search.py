import csv
import numpy as np
import pandas as pd
import os
from Bio import SeqIO
import gzip, pickle
import random
import shutil
from tqdm import tqdm
tqdm.monitor_interval = 0
from ete2 import NCBITaxa
import itertools
import subprocess
import glob
from time import time
from pprint import pprint
from scipy import stats

PATH_TO_GENOMES = 'input'
PARENT_LEVEL_OF = {"order": "class",
                   "family": "order",
                   "genus": "family",
                   "species": "genus"}

def get_desired_ranks(taxid, desired_ranks):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())

    return [ranks2lineage.get(rank, '0') for rank in desired_ranks]

def get_parent2genome(parent_level="genus"):

    parent2genome = dict()

    for fasta in os.listdir(os.path.join(PATH_TO_GENOMES)):
        if fasta.endswith('.txt'):
            continue
        key = get_desired_ranks(fasta.split('.')[0], [parent_level])[0]
        if not parent2genome.has_key(key):
            parent2genome[key] = []
        parent2genome[key].append(fasta)
    return parent2genome

def parse_square_matrix(df):
    # Get 75th percentile of current matrix
    # This will be used as a threshold for meta-jaccard estimation
    series = pd.Series(df.values.ravel())
    summary_stats = series.describe()

    # Threshold inbetween 50th and 75th percentile
    threshold = (summary_stats[5] + summary_stats[6]) / 2

    # Get Jaccard distances for every genome pair in the matrix
    # 0 = Complete similarity; 1 = Completely dissimilarity
    # Output dictionary - (Genome 1, Genome 2): Jaccard Distance
    rows_as_pandas_series = []
    column_names = []
    for index, row in df.iterrows():
        if index == 0:
            column_names = row.index
        rows_as_pandas_series.append(row)
    
    col_index_to_name = dict()
    for index, name in enumerate(column_names):
        col_index_to_name[index] = name
    
    jaccard_distances = dict()
    for index, row in enumerate(rows_as_pandas_series):
        for index_col, distance in enumerate(row):
            if index_col != 0:
                jaccard_distances[ (int(row[0]), int(col_index_to_name[index_col])) ] = distance
    
    return (jaccard_distances, threshold)

def get_all_jaccard_matrices(path_to_matrices):
    start_time = time()

    jaccard_matrices_and_threshold_pairs = dict()
    for kmer_matrix in tqdm(os.listdir(path_to_matrices)):
        kmer_index = int((kmer_matrix.split('_')[-1]).split('.')[0])

        kmer_matrix = os.path.join(path_to_matrices, kmer_matrix)
        df = pd.read_csv(kmer_matrix, sep=";", index_col=1)
        (jaccard_distance, threshold) = parse_square_matrix(df)

        jaccard_matrices_and_threshold_pairs[kmer_index] = (jaccard_distance, threshold)

    end_time = time()
    print("All Jaccard Matrices loaded in: {} seconds").format(end_time - start_time)

    return jaccard_matrices_and_threshold_pairs

def search_for_best_k(parent2child, key1, key2, kmer_values,
                      jaccard_matrices_and_threshold_pairs):

    for kmer_len in kmer_values:

        (jaccard_distances, meta_jaccard_threshold) = jaccard_matrices_and_threshold_pairs[kmer_len]

        rows = parent2child[key1]
        columns = parent2child[key2]
        total_sum = 0

        for i, row_cell in enumerate(rows):
            for j, col_cell in enumerate(columns):
                row_cell_name = int(row_cell.split('.')[0])
                col_cell_name = int(col_cell.split('.')[0])
                try:
                    total_sum += jaccard_distances[ (row_cell_name, col_cell_name) ]
                except:
                    total_sum += jaccard_distances[ (col_cell_name, row_cell_name) ]
        
        meta_jaccard_value = total_sum / (len(rows) * len(columns))

        # If meta jaccard value is greater than 75th percentile dissimilarity of the
        # current full k-mer matrix, then the kmer value is a good separating value between
        # the pair of parent taxon nodes being compared
        if (meta_jaccard_value > meta_jaccard_threshold):
            return kmer_len

    return kmer_values[-1]


if __name__ == "__main__":

    # Load jaccard distance matrices for each kmer length
    path_to_matrices = os.path.join("..", "jaccard_distances_calculated")

    jaccard_matrices_and_threshold_pairs = get_all_jaccard_matrices(path_to_matrices)
    kmer_values = (jaccard_matrices_and_threshold_pairs.keys())
    kmer_values.sort()

    parent_levels = ["order", "family", "genus", "species"]

    kmer_characterization = dict()

    for parent_level in parent_levels:

        parent2child = get_parent2genome(parent_level)
        permutations = list(itertools.combinations(parent2child.keys(), 2))

        for vals in permutations:
            common_ancestor_0 = get_desired_ranks(vals[0], [PARENT_LEVEL_OF[parent_level]])[0]
            common_ancestor_1 = get_desired_ranks(vals[1], [PARENT_LEVEL_OF[parent_level]])[0]

            print "Common ancestor of {} and {} is {} - Level {}".format(vals[0], vals[1], common_ancestor_0,
                                                                         PARENT_LEVEL_OF[parent_level])

            if common_ancestor_0 == common_ancestor_1:
                if not kmer_characterization.has_key(common_ancestor_0):
                    kmer_characterization[common_ancestor_0] = []
            
                kmer_len = search_for_best_k(parent2child, vals[0], vals[1], kmer_values,
                                             jaccard_matrices_and_threshold_pairs)
                kmer_characterization[common_ancestor_0].append(kmer_len)
        
        print "-----------------------------------------------------"
    
    for k, v in kmer_characterization.items():
        print"{}: {}".format(k, sum(v) / len(v))
    
