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

	hgeXY_dist = genetree_copy.get_distance("hge_x", "hge_y")
	hliX_dist = genetree_copy.get_distance("hli", "hge_x")
	hliY_dist = genetree_copy.get_distance("hli", "hge_y")
	

	if hgeXY_dist < hliX_dist and hgeXY_dist < hliY_dist:
		return "hgeXY"
	elif hliX_dist < hgeXY_dist and hliX_dist < hliY_dist:
		return "hliX"
	elif hliY_dist < hliX_dist and hliY_dist < hgeXY_dist:
		return "hliY"

def calc_topology_asymmetry(genetrees):

	hgeXY_count, hliX_count, hliY_count = 0, 0, 0

	for tree in genetrees:
	
		topology = get_tree_topology(tree)

		if topology == "hgeXY":
			hgeXY_count += 1
		elif topology == "hliX":
			hliX_count += 1
		elif topology == "hliY":
			hliY_count += 1

	return hgeXY_count, hliX_count, hliY_count 

hgeXY, hliX, hliY = calc_topology_asymmetry(trees)
print(hgeXY, hliX, hliY)


		
