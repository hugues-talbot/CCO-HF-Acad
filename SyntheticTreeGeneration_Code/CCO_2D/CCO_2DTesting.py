# -*- coding: utf-8 -*-
"""
Created on Fri May 06 14:28:15 2016

@author: cjaquet
"""

##HT not sure what this is for
##HT 

import numpy as np
import matplotlib.pyplot as plt
import Kamiya as kami
#import NodeClass as nclass
import copy
from pylab import figure, gca, Line2D
import CCO_2D as cco

fd = open('./Results/CCO_beta2_output.txt','w') # open the result file in write mode
old_stdout = sys.stdout   # store the default system handler to be able to restore it    
sys.stdout = fd # Now your file is used by print as destination 

first_tree = cco.cco_function(5, "NodeClass")
sec_tree = cco.cco_function(5, "NodeClassBis")

first_tree.printing_full()
sec_tree.printing_full()

sys.stdout=old_stdout # here we restore the default behavior
fd.close() # to not forget to close your file

if False:
    Q_perf = 8.33e3
    N_term = 30.		
    Q_term = Q_perf / N_term
    P_drop = 1.33e7 - 8.38e6
    viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm)
    
    # tree 
    tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity)
    
    #source point 		
    area_center = np.array([100.,50.])
    area_radius = 50
    area_descptr = [area_center, area_radius]
    area = np.pi*area_radius**2
    root_position = np.array([100,100])
    root_node = nclass.Node(0,root_position, Q_perf, -1)
    root_node.set_child_0_index(1)
    tree.add_node(root_node)
    
    
    #first segment end
    
    first_node_position = np.array([100,60]) #first_segmt_end(area)# question: should it be over a threshold ?
    print "first segment end point", first_node_position
    first_node = nclass.Node(1, first_node_position, Q_perf,0)
    first_node.set_child_0_index(3)
    first_node.set_child_1_index(2)
    ##first_node.set_betas(np.array([0.6,0.8]))
    tree.add_node(first_node)
            
    second_node_pos = np.array([80,40])
    second_node = nclass.Node(2,second_node_pos, Q_perf,1)
    second_node.set_child_0_index(4)
    second_node.set_child_1_index(5)
    tree.add_node(second_node)
    #
    second_bis_pos = np.array([120,40])
    sc_bis = nclass.Node(3,second_bis_pos, Q_perf, 1)
    tree.add_node(sc_bis)
    #
    third_node_pos = np.array([60,20])
    third_node = nclass.Node(4, third_node_pos, Q_perf, 2)
    tree.add_node(third_node)
    #
    fourth_node_pos = np.array([90,20])
    fourth = nclass.Node(5, fourth_node_pos, Q_perf, 2)
    tree.add_node(fourth)
    
    tree.update_flow()
        
    cco.plot_tree(tree, area_descptr)
    
    #tree.balancing_ratios(5)
    
    
    #cco.plot_tree(tree, area_descptr)
    
    #####testing intersections###################
    store_locations = []
    store_good_locations = []
    for i in range (50):
        print i
        location_1 = cco.first_segmt_end(area, area_descptr)#np.array([126.,33.])#
        location_2 = cco.first_segmt_end(area, area_descptr)#np.array([102.,65.]) #
        k_term = tree.get_k_term()
        #for j in range (k_term):
        if (tree.check_intersection(0, location_1, location_2, 0.96)):
            print "locations", location_1, location_2
            #init = nclass.Node(k_term + 1, location_1,Q_perf, -1)
            #tree.add_node(init)
            #finish = nclass.Node(k_term+1, location_2,Q_perf, k_term)
            #tree.add_node(finish)
            store_good_locations.append([location_1, location_2])
        else:
            store_locations.append([location_1, location_2])
            
            
    def plot_tree_and_list(tree, area_descptr, list_a, list_b):
        fig = plt.figure(figsize=(8,8))#figsize=(8,8)
        ax = fig.add_subplot(111)
        ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[1], color = 'r', alpha = 0.5))
        for sgmt in tree.nodes:
            if (sgmt.parent() >= 0):
                distal = sgmt.coord
                proximal = tree.get_node(sgmt.parent()).coord
                radius = 0.96#*10#*10 #tree.get_radius(sgmt.index) *10
                #ax.plot([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'b')
                ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'b'))
                
        ax.set_xlim(area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5)
        ax.set_ylim(area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5) 
        for i in list_a:
            #plt.plot([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth = 0.96, color = 'b')
            ax.add_line(Line2D([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth=0.96, color = 'b'))
        for i in list_b:
            #plt.plot([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth = 0.96, color = 'r')
            ax.add_line(Line2D([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth=0.96, color = 'r'))
        fig.savefig("./intersections_test.png")
        fig.show()        
    
    plot_tree_and_list(tree, area_descptr, store_good_locations, store_locations)    
    #cco.plot_tree(tree, area_descptr)
