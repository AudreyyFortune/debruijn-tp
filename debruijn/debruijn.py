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

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

def read_fastq(file):
	"""
	Read fastq file.
	input : fastq
	output : generateur de sequence
	"""
	with open(file) as fichier:
		for line in enumerate(fichier):
			yield(next(fichier).strip()) # on saute la ligne commençant par "@".
			next(fichier) # on saute la ligne avec les 'JJJJ'.
			next(fichier) # on saute la ligne commençant par "@".


def cut_kmer(seq, taille):
	"""
	Cut a sequence in k-mers.
	input : sequence, k-size
	output : generateur de k-mers
	"""
	for n in range(len(seq) - taille + 1):  
		yield(seq[n:n+taille])


def build_kmer_dict(file, taille):
	"""
	"""
	dico = {}
	for seq in read_fastq(file):
		for kmer in cut_kmer(seq, taille):
			if kmer not in dico:
				dico[kmer] = 0
			dico[kmer] += 1
	print(dico)
	return(dico)



if __name__ == '__main__':
	main()
	#for i in read_fastq("../data/eva71_two_reads.fq"):
	#	print(i)
	#	for j in cut_kmer(i,4):
	#		print(j)
	#print(build_kmer_dict("../data/eva71_two_reads.fq",4))
	args = get_arguments()
	build_kmer_dict(args.fastq_file, args.kmer_size)