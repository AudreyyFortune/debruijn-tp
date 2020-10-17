#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
from graphviz import Digraph
import matplotlib.pyplot as plt

__author__ = "Fortune Audrey"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Fortune Audrey"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Fortune Audrey"
__email__ = "audefortune@free.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


#==============================================================
# Part 1
#==============================================================

def read_fastq(file):
	"""
	Read fastq file.
	Parameters:
		input : fastq (name of the file)
		return : generator of sequence
	"""
	with open(file) as fichier:
		for line in enumerate(fichier):
			yield(next(fichier).strip()) # on saute la ligne commençant par "@".
			next(fichier) # on saute la ligne avec les 'JJJJ'.
			next(fichier) # on saute la ligne commençant par "@".


def cut_kmer(seq, taille):
	"""
	Cut a sequence in k-mers.
	Parameters:
		input : - seq (a sequence to be cut into k-mers)
				- taille (k-size, size of the k-mers)
		return : generator of k-mers
	"""
	for n in range(len(seq) - taille + 1):  
		yield(seq[n:n+taille])


def build_kmer_dict(file, taille):
	"""
	Build a k-mers dictionnary.
	Parameters:
		input : - fastq (name of the file)
				- taille (k-size)
		return : dictionnary (keys : taille, values : occurancy)
	"""
	dico = {}
	for seq in read_fastq(file): # seq
		for kmer in cut_kmer(seq, taille): # k-mers
			if kmer not in dico:
				dico[kmer] = 0
			dico[kmer] += 1
	return(dico)


def build_graph(dico):
	"""
	Creation of graph of all k-mer prefixes and suffixes.
	Parameters:
		input : dictionnary of kmers and their occurencies
		return : graph of kmers with links between prefix and suffix
	"""
	G = nx.DiGraph()
	for k in dico:
		prefix = k[:-1]
		suffix = k[1:]
		G.add_node(prefix) # add_node -> ajoute un noeud objet au graphe. 
		G.add_node(suffix)
		G.add_edge(prefix, suffix, weight = dico[k]) # add_edge -> ajoute un lien entre prefix et suffix. 
	#nx.draw(G,with_labels = True) #pour tracer
	#plt.show()
	return G


#==============================================================
# Part 2
#==============================================================

def get_starting_nodes(G):
	"""
	Gives entry nodes of the graph.
	Parameters:
		input : graph
		output : list of all entry nodes
	"""
	starting_node = []
	for node in G.nodes:
		if not G.in_edges(node): # in_edges -> renvoie une liste des bords entrants.
			starting_node.append(node)
	return starting_node


def get_sink_nodes(G):
	"""
	Gives output nodes of the graph.
	Parameters:
		input : graph
		output : list of all exit nodes
	"""
	sink_node = []
	for node in G.nodes:
		if not G.out_edges(node): # out_edges -> renvoie une liste des bords sortants.
			sink_node.append(node)
	return sink_node


def get_contigs(G, starting_node, sink_node):
	"""
	Determine each path of the graph between a starting node and a strating node, named "contig".
	Parameters:
		input: - graph
			   - starting_node(list of all entry nodes)
			   - sink_node(list of all exit nodes)
		output: list of tupple containing contig and contig size
	"""
	contigs = []
	for start_n in starting_node:
		for exit_n in sink_node:
			for chemin in nx.all_simple_paths(G, start_n, exit_n): # all_simple_paths -> genere tous les chemins simple dans le graphe.
				cont = chemin[0]
				for i in range(1,len(chemin)):
					cont += chemin[i][-1]
				contigs.append((cont,len(cont)))
	return contigs
   

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
     

def save_contigs(contigs, output):
	"""
	Saves and writes an output file containing contigs in fasta format.
	Parameters:
		input: - contigs(list of tuple(contig, contig size))
			   - output (name of the output file)
	"""
	with open(output, "w") as filout:
		for index, element in enumerate(contigs):
			filout.write(">contig_" + str(index) + " len=" + str(element[1]) + "\n")
			filout.write(fill(element[0]) + "\n")


    

#==============================================================
# Part 3
#==============================================================

def std():
	"""
	"""
	pass


def path_average_weight():
	"""
	"""
	pass
	
	
def remove_paths():
	"""
	"""
	pass
	
def select_best_path():
	"""
	"""
	pass
	
def solve_bubble():
	"""
	"""
	pass

def simplify_bubbles():
	"""
	"""
	pass



#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # K-mers dictionnary
    dico_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    # Graph de Bruijn 
    graphe = build_graph(dico_kmer)
    # Starting nodes of graph
    start_node = get_starting_nodes(graphe)
    # Sinking nodes of graph
    exit_node = get_sink_nodes(graphe)
    # Get contigs
    contigs = get_contigs(graphe, start_node, exit_node)
    # Save contigs
    save_contigs(contigs, args.output_file)
    
    
if __name__ == '__main__':
	main()