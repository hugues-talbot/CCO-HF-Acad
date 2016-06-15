# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 13:48:17 2016

@author: cjaquet
"""

##HT It seems this code is only for testing, please confirm

import numpy as np
import matplotlib.pyplot as plt
import Kamiya as kami
import NodeClass as nclass
import CCO_2D as cco

Q_perf = 8.33e3
N_term = 250.		
Q_term = Q_perf / N_term
P_drop = 1.33e7 - 8.38e6
viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm)

######## Testing Kamiya's algorithm #################

##initialization

# tree 
tree = nclass.Tree([], 24, Q_perf, P_drop, viscosity)

#source point 		
area_center = np.array([100.,50.])
area_radius = 50
area_descptr = [area_center, area_radius]
area = np.pi*area_radius**2
root_position = np.array([95,75])#100,100
root_node = nclass.Node(0,root_position, Q_perf, -1)
tree.add_node(root_node)


#first segment end
#generate random position > dist_crit

first_node_position = np.array([130,25])#100,20#first_segmt_end(area)# question: should it be over the d_tresh?
print "first segment end point", first_node_position
first_node = nclass.Node(1, first_node_position, Q_perf,0)
#first_node.set_betas(np.array([0.6,0.8]))
tree.add_node(first_node)


test = np.array([140,80])#75,40# np.array([75,40])#get_new_location(tree, area_descptr, N_term)
print "new location", test
#second_node = nclass.Node(2, test, Q_perf, 1)
#tree.add_node(second_node)
radius_ori = tree.get_radius(1)
r = np.ones(3)*radius_ori
print "initial radii", r
f = np.array([5666,2833,2833])#Q_perf, Q_perf*1./4., Q_perf*1./4.
print "initial flow", f
c0 = root_position
c1 = first_node_position
c2 = test
print "initial positions", c0,c1,c2
iter_max = 100
tolerance = 0.01
convergence, res = kami.kamiya_single_bif_opt(r, f, c0, c1, c2, iter_max, tolerance, True)
gamma = 3.
ind = 0
trees = []
for i in res:
    if ind == 0 or ind == len(res)-1: 
        tre = nclass.Tree([], 24, Q_perf, P_drop, viscosity)
        
        tre.add_node(root_node)
        
        bif_node = nclass.Node(1, np.array([i[0], i[1]]), Q_perf,0)
        betas = nclass.calculate_betas(i[2][1]/i[2][2],gamma)
        bif_node.set_betas(betas)
        bif_node.set_old_child_index(2)
        bif_node.set_young_child_index(3)
        tre.add_node(bif_node)
        
        old_child = nclass.Node(2, first_node_position, Q_perf/2.,1)
        tre.add_node(old_child)
        
        new_child = nclass.Node(3, test, Q_perf/2.,1)
        tre.add_node(new_child)
        trees.append(tre)
        #cco.plot_tree(tre, area_descptr, str(ind))
    ind = ind + 1
    
cco.plot_trees(trees, area_descptr)
