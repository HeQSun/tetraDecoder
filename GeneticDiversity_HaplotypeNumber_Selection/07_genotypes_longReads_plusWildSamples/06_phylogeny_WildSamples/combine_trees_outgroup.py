#!/usr/bin/env python

"""                                                                                 
python merge_seqs.py list_samples.txt 'alig_*'
"""                                                                                 

from ete3 import Tree
from Bio import SeqIO
from Bio.Seq import Seq

import glob
import sys

# list of trees
listTrees = [line.strip() for line in open(sys.argv[1])]

newtrees = []
for line in listTrees:
	# line_tree=open(line, 'r')
	tree_line = Tree(line)
	# set outgroup:
	tree_line.set_outgroup(str(sys.argv[3]))
	newtrees += [tree_line.write()]
	# print(newtrees)

# Print trees in a new file
rootedtrees = open(sys.argv[2], 'w')
for tree in newtrees:
	rootedtrees.write(tree + "\n")


