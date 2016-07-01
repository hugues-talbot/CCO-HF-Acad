# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import matplotlib.pyplot as plt
import Kamiya as kami
import NodeClass as nclass
import CCO_2DFunctions as cco_2df
import copy
from pylab import figure, gca, Line2D
import random
import os
import cPickle as pickle
import sys  # Need to have acces to sys.stdout
import time
from multiprocessing import Pool, TimeoutError




############# Visualisation tools ####################

def plot_tree(tree, area_descptr, name):
    fig = plt.figure(figsize=(8,8))#figsize=(8,8)
    ax = fig.add_subplot(111)
    #colors=['k', 'b', 'g', 'm', 'r']
    # labels for cet ['r','m','g','b']
    ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[1], color = 'r', alpha = 0.5))
    all_radii = []
    leaves_radii = []
    inv_length_fac = 1. / tree.length_factor
    for sgmt in tree.nodes:
        if (sgmt.parent() >= 0):
            distal = sgmt.coord
            proximal = tree.get_node(sgmt.parent()).coord
            radius = tree.get_radius(sgmt.index) *inv_length_fac #*10
            all_radii.append(radius)
            if sgmt.is_leaf():
                leaves_radii.append(radius)               
            #ax.plot([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'b')
            ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'k'))#colors[sgmt.label]
      
    ax.set_xlim(area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5)
    ax.set_ylim(area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5)
    outfile_path1 = "./"+name+"all.txt"
    outfile1 = open(outfile_path1,"w")
    for i in all_radii:
        outfile1.write(str(i) + "\n")
        
    outfile_path2 = "./"+name+"leaves.txt"
    outfile2 = open(outfile_path2,"w")
    for i in leaves_radii:
        outfile2.write(str(i) + "\n")
    
    fig.savefig(name+".png")
    fig.savefig(name+".pdf")
    fig.show()

def plot_trees(trees, area_descptr):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    colors = ['b', 'r']
    ind = 0
    for tree in trees:      
        for sgmt in tree.nodes:
            if (sgmt.parent() >= 0):
                distal = sgmt.coord
                proximal = tree.get_node(sgmt.parent()).coord
                radius = tree.get_radius(sgmt.index) *10
                
                ax.plot([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = colors[ind], alpha= 0.4)
        ind = ind +1   
    ax.set_xlim(0,200)#area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5
    ax.set_ylim(0,200)#area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5

    return fig
    
##################################################################################
############# Karch algo : CCO ####################

timing = True
store_data = False

if timing:
    debut = time.time()
if store_data:
    fd = open('./Results/CCO_test.txt','w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it    
    sys.stdout = fd # Now your file is used by print as destination 
    

if True:

    NTerm = 4000
    seed = 42

    np.random.seed(seed)
    process_nb = 8 
 
    #### Parameters to define: ##
    ## About tree
      
    Q_perf = 8.33e3
    N_term = NTerm
    Q_term = Q_perf / N_term
    P_drop = 1.33e7 - 7.98e6 # when Nterm = 4 000, the P_drop is 1.33e7 -7.98e6 #when =Nterm=250 :1.33e7 - 8.38e6
    viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm) #3.6 cp =3.6e-3 Pa.s = 3.6e-9 MPa.s 


    N_con = 20
    N_con_max = 40
    
    # About  convexe perfusion surface : defines a disc surface 
    area_center = np.array([100.,50.])
    area_radius = 50
    area_descptr = [area_center, area_radius]
    area = np.pi*area_radius**2    
    
    
    #### initialization ##    
    # tree
    store_cet = []
    tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity, area, area_radius)
     
    # source point : define the root position
    root_position = np.array([100,100])
    root_node = nclass.Node(0,root_position, Q_perf, -1)
    root_node.set_child_0_index(1)
    root_node.set_label(0)
    tree.add_node(root_node)    
    
    #first segment end: randomly picked inside perfusion surface
    tree.update_length_factor()
    first_node_position = cco_2df.first_segmt_end(area, area_descptr)
    first_node = nclass.Node(1, first_node_position, Q_perf,0)
    first_node.set_label(0)
    tree.add_node(first_node)

    tree.update_flow()
    tree.update_length_factor()
    tree.depthfirst_resistances(0)

            
    process_nb = 4            
    while tree.get_k_term() < N_term: 

        success, new_child_location, d_tresh = cco_2df.get_new_location(tree, area_descptr, N_term)
        if (success == False):
            #print "impossible to satisfy distance criteria", "d_tresh", d_tresh
            ## then try again with
            break

        cet = []#[] for i in range (N_con_max)c=

        adding_location = False
        added_location = []
        
        dtype_r=[("convgce", int), ("volume", float), ("betas_and_branching_location", float, (2, 2)),("old_child_index", int)]
                
        test_N_con_max = False
        
        # test closest neighbors
        neighbors = tree.find_neighbors(new_child_location, N_con)

        args = [[tree, neighbors[i],new_child_location] for i in range (len(neighbors))]

        #process all neighbors connection test batch by batch
        while (len(cet) < N_con) and (len(cet) < len(neighbors)):           
            end_lim = len(cet) + process_nb if (len(cet) + process_nb < len(neighbors)) else len(neighbors)            
            pool = Pool(processes =  process_nb)               
            res = pool.map(cco_2df.test_connection_list,args[len(cet): end_lim]) 
            cet = cet + res
            pool.close()

        cet_filtered = filter(None,cet)
        cet_values = np.array(cet_filtered, dtype_r)
        
        #if there is a candidate that converges
        if (np.sum(cet_values['convgce']) > 1) or (np.sum(cet_values['convgce']) > 0 and tree.get_k_term() == 1):        
            cet_sel = cet_values[cet_values['convgce']>0]
            cet_sorted = np.sort(cet_sel, order = "volume")
            cet_final=cet_sorted[0]
            adding_location = True
            added_location.append(cet_final.tolist()[1:])

        # test extra neighbors if no connection candidate has fullfilled constraints
        else: 
	    if len(neighbors) == N_con:  # if there are at least N_con neighbors there might be extra neighbors to test with
	      #print "testing extra"
	      test_N_con_max = False
	      neighbors = tree.find_neighbors(new_child_location, N_con_max)
	      extra_neighbs = neighbors[N_con:N_con_max]
	      argsb = [[tree, extra_neighbs[i],new_child_location] for i in range (len(extra_neighbs))]
	      #process all neighbors connection test batch by batch
	      while (len(cet) < N_con_max) and (len(cet) < N_con + len(extra_neighbs)):
               end_lim = (len(cet) + process_nb - N_con) if (len(cet) + process_nb - N_con < len(extra_neighbs)) else len(extra_neighbs)  
               poolb = Pool(processes =  process_nb)               
               res = poolb.map(cco_2df.test_connection_list,argsb[len(cet) - N_con: end_lim]) 
               cet = cet + res 
               poolb.close()

	      cet_filtered = filter(None,cet)
	      cet_values = np.array(cet_filtered, dtype_r)   
	      if (np.sum(cet_values['convgce']) > 1) or (np.sum(cet_values['convgce']) > 0 and tree.get_k_term() == 1):
               cet_values = np.array(cet_filtered, dtype_r)
               cet_sel = cet_values[cet_values['convgce']>0]
               cet_sorted = np.sort(cet_sel, order = "volume")
               cet_final=cet_sorted[0]
               adding_location = True
               added_location.append(cet_final.tolist()[1:])

        

        if (adding_location): # optimal connection found!

            store_cet.append(filter(None,cet))

            opt = added_location[-1]
            ante_tree = copy.deepcopy(tree)
            
            if (tree.add_connection(opt[2], new_child_location, opt[1][1], opt[1][0])):
                branching_added_node = tree.nodes[len(tree.nodes) - 2]       
                added_node = tree.nodes[len(tree.nodes) - 1]
                label = int(len(tree.nodes) / 100.)             
                if (label == 5):
                    label = 4
                added_node.set_label(label)
                branching_added_node.set_label(label)

                print "connection added on tree up node", opt[2]
                print "k term is now ", tree.get_k_term()

                last_tree = copy.deepcopy(tree)
                if tree.get_k_term() % 50  == 0:
                    plot_tree(tree, area_descptr, "./Results/tree_Nt%i_s%i_final" %(tree.get_k_term(),seed))#tree_stored[-1]
            else:
                print "failed to add connection on tree"
        else:
            print "location doesn't provide an optimal connection, testing newwwwwwwww location"

        #keep going until reach Nterm!

    plot_tree(tree, area_descptr, "./Results/tree_Nt%i_s%i_final" %(tree.get_k_term(),seed))#tree_stored[-1]





if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"
    

