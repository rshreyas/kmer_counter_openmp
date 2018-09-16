Author: Shreyas Ramesh

This folder contains all files necessary to run the kmer counting procedure mentioned in the Report "Parallel Multi-set K-mer Counting on OpenMP for Hierarchical Taxonomy Classification" by Shreyas Ramesh.


DIRECTORIES:
	fasta_reader - fasta file parser used by all downstream k-mer counting strategies

	input - Contains all 210 genome FASTA files used for k-mer counting

	jaccard_distances_calculated - Jaccard Distances between every genome pair at k-mer sizes ranging from 8 to 31.

	graph_generation - Jupyter Notebook that contains code that was used to generate graphs in the report
	
	optimal_kmer_search - Contains a script (optimal_kmer_size_search.py) to find optimal k-mer size from Jaccard abundance matrices.



FILES:
	execute.sh - Cleans, recompiles and executes the kmer_counter program

	kmer_counter.cpp - Main logic for k-mer counting, along with the OpenMP pragmas

	Makefile - Compiles the kmer_counting program, links the fasta_reader object file and uses the correct flags

	optimal_kmer_size_search.py - Python implementation that finds the best value of 'k' for the dataset.
				      This program uses the Jaccard distances within the simka_results folder.
				      Run with command: python optimal_kmer_size_search.py

	Graphs for Parallel Processing Assignment.ipynb - Jupyter Notebook that contains logic to generate all graphs used
							  within the report.


REQUIREMENTS: GNU g++ compiler with OpenMP 3.0 and above.
	      This program was tested on ironclaw1-4 and a linux operating system.


Compile by running "execute.sh"


All parameters in this project can be changed in the "driver" procedure within kmer_counter.cpp


