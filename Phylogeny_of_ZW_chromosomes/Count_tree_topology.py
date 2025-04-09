######## this script is modified from https://github.com/mhibbins/RumexComparative/blob/main/buceph_origins/count_sexlinked_trees.py

import sys
from ete3 import Tree
import numpy as np
import random
import copy

trees = []

with open(sys.argv[1], 'r') as treefile:
	for line in treefile:
		trees.append(Tree(line))

def get_tree_topology(genetree):

	genetree_copy = copy.deepcopy(genetree)

	hsaZW_dist = genetree_copy.get_distance("Hsa", "Hsab")
	hinZ_dist = genetree_copy.get_distance("Hni", "Hsa")
	hinW_dist = genetree_copy.get_distance("Hni", "Hsab")
	

	if hsaZW_dist < hinZ_dist and hsaZW_dist < hinW_dist:
		return "hsaZW"
	elif hinZ_dist < hsaZW_dist and hinZ_dist < hinW_dist:
		return "hinZ"
	elif hinW_dist < hinZ_dist and hinW_dist < hsaZW_dist:
		return "hinW"

def calc_topology_asymmetry(genetrees):

	hsaZW_count, hinZ_count, hinW_count = 0, 0, 0

	for tree in genetrees:
	
		topology = get_tree_topology(tree)

		if topology == "hsaZW":
			hsaZW_count += 1
		elif topology == "hinZ":
			hinZ_count += 1
		elif topology == "hinW":
			hinW_count += 1

	return hsaZW_count, hinZ_count, hinW_count 

hsaZW, hinZ, hinW = calc_topology_asymmetry(trees)
print(hsaZW, hinZ, hinW)


		
