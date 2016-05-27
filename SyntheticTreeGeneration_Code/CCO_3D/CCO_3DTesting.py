# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:28:37 2016

@author: cjaquet
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, gca, Line2D
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
import copy

import Kamiya_3D as kami
import NodeClass3D as nclass
#import CCO_3D as cco
import CCO_3DFunctions as cco_3df




Q_perf = 8.33e3
N_term = 30.		
Q_term = Q_perf / N_term
P_drop = 1.33e7 - 8.38e6
viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm)

# tree 
tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity)

#source point 		
area_center = np.array([100.,50.,0])
area_radius = 50
area_descptr = [area_center, area_radius]
area = np.pi*area_radius**2

volume_descptr = [area_center, area_radius]
root_position = np.array([100,100,0])
root_node = nclass.Node(0,root_position, Q_perf, -1)
root_node.set_child_0_index(1)
tree.add_node(root_node)


#first segment end

first_node_position = np.array([100,60,0]) #first_segmt_end(area)# question: should it be over a threshold ?
print "first segment end point", first_node_position
first_node = nclass.Node(1, first_node_position, Q_perf,0)
first_node.set_child_0_index(3)
first_node.set_child_1_index(2)
##first_node.set_betas(np.array([0.6,0.8]))
tree.add_node(first_node)
        
second_node_pos = np.array([80,40,0])
second_node = nclass.Node(2,second_node_pos, Q_perf,1)
second_node.set_child_0_index(4)
second_node.set_child_1_index(5)
tree.add_node(second_node)
#
second_bis_pos = np.array([120,40,0])
sc_bis = nclass.Node(3,second_bis_pos, Q_perf, 1)
tree.add_node(sc_bis)
#
third_node_pos = np.array([60,20,0])
third_node = nclass.Node(4, third_node_pos, Q_perf, 2)
tree.add_node(third_node)
#
fourth_node_pos = np.array([90,20,0])
fourth = nclass.Node(5, fourth_node_pos, Q_perf, 2)
tree.add_node(fourth)

tree.update_flow()

#####testing intersections###################
store_locations = []
store_good_locations = []
for i in range (1):
    print i
    location_1 = cco_3df.first_segmt_end(volume_descptr)#np.array([58.081,20,0])
    location_2 = cco_3df.first_segmt_end(volume_descptr)#np.array([80,60,20])
    location_1[2]=0.
    #location_2[2]=0.
    k_term = tree.get_k_term()
    #for j in range (k_term):
    radii = np.ones(3)*0.96
    if (tree.check_intersection(0, location_1, location_2, radii)):
        print "locations", location_1, location_2
        #init = nclass.Node(k_term + 1, location_1,Q_perf, -1)
        #tree.add_node(init)
        #finish = nclass.Node(k_term+1, location_2,Q_perf, k_term)
        #tree.add_node(finish)
        store_good_locations.append([location_1, location_2])
    else:
        store_locations.append([location_1, location_2])
        
        
def plot_tree_and_list(tree, area_descptr, list_a, list_b):
    wid=8
    hei=8
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    #ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[1], color = 'r', alpha = 0.5))
    for sgmt in tree.nodes:
        if (sgmt.parent() >= 0):
            distal = sgmt.coord
            proximal = tree.get_node(sgmt.parent()).coord
            radius = 0.96#*10#*10 #tree.get_radius(sgmt.index) *10
            #ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'b'))
            verts = [zip([distal[0], proximal[0]],[distal[1], proximal[1]],[distal[2], proximal[2]])]
            tri = a3.art3d.Poly3DCollection(verts)
            tri.set_color('b')
            tri.set_linewidth(radius*72*hei/(area_descptr[1]*3.))
            ax.add_collection3d(tri)
                
    ax.set_xlim(area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5)
    ax.set_ylim(area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5) 
    for i in list_a:
        #ax.add_line(Line2D([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth=0.96, color = 'b'))
        verts = [zip([i[0][0], i[1][0]],[i[0][1], i[1][1]],[i[0][2], i[1][2]])]
        tri = a3.art3d.Poly3DCollection(verts)
        tri.set_color('g')
        tri.set_alpha(0.5)
        tri.set_linewidth(radius*72*hei/(area_descptr[1]*3.))
        ax.add_collection3d(tri)
        
    for i in list_b:
        #ax.add_line(Line2D([i[0][0], i[1][0]],[i[0][1], i[1][1]], linewidth=0.96, color = 'r'))
        verts = [zip([i[0][0], i[1][0]],[i[0][1], i[1][1]],[i[0][2], i[1][2]])]
        tri = a3.art3d.Poly3DCollection(verts)
        tri.set_color('r')
        tri.set_alpha(0.5)
        tri.set_linewidth(radius*72*hei/(area_descptr[1]*3.))
        ax.add_collection3d(tri)
    
    
    ax.set_xlim(50,100)#area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5
    ax.set_ylim(50,100)#area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5
    ax.set_zlim(-50,50)#area_descptr[0][2] - area_descptr[1]*1.5,area_descptr[0][2]+ area_descptr[1]*1.5
    plt.savefig("./Results/intersections_test.png")
    plt.show()
    #fig.show()        
    
plot_tree_and_list(tree, volume_descptr, store_good_locations, store_locations)