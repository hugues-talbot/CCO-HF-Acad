# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:10:48 2016

@author: jaquetc
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import matplotlib.pyplot as plt
import Kamiya_3D as kami
import NodeClass3D as nclass
import CCO_3DFunctions as cco_3df
import copy

import random
import os
import cPickle as pickle
import sys  # Need to have acces to sys.stdout
import time

from pylab import figure, gca, Line2D
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d
import pylab as pl

##########################################
def plot_tree(tree, vol_descptr, name, factor_radius):
    wid=16
    hei=16
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    #colors=['k', 'b', 'g', 'm', 'r']
    # labels for cet ['r','m','g','b']
    #ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[1], color = 'r', alpha = 0.5))
    #p = Circle((vol_descptr[0][0], vol_descptr[0][1]), vol_descptr[1], color = 'r', alpha = 0.5)
    #ax.add_patch(p)
    #art3d.pathpatch_2d_to_3d(p, z=vol_descptr[0][2], zdir="z")
    
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = vol_descptr[0][0]+ vol_descptr[1] * np.outer(np.cos(u), np.sin(v))
    y = vol_descptr[0][1]+ vol_descptr[1] * np.outer(np.sin(u), np.sin(v))
    z = vol_descptr[0][2]+ vol_descptr[1] * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1) #,
    
    for sgmt in tree.nodes:
        if (sgmt.parent() >= 0):
            distal = sgmt.coord
            proximal = tree.get_node(sgmt.parent()).coord
            radius = tree.get_radius(sgmt.index) #* factor_radius 
            #ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius, color = 'b'))#colors[sgmt.label]
            verts = [zip([distal[0], proximal[0]],[distal[1], proximal[1]],[distal[2], proximal[2]])]
            tri = a3.art3d.Poly3DCollection(verts)
            tri.set_color('k')
            tri.set_linewidth(radius*72*hei/(vol_descptr[1]*3.) )
            ax.add_collection3d(tri)
            
    bound = vol_descptr[1]*1.5
    ax.set_xlim(vol_descptr[0][0] - bound,vol_descptr[0][0]+ bound)
    ax.set_ylim(vol_descptr[0][1] - bound,vol_descptr[0][1]+ bound)
    ax.set_zlim(vol_descptr[0][2] - bound,vol_descptr[0][2]+ bound)
    ax.grid(False)
    plt.savefig(name+".png")
    plt.savefig(name+".pdf")
    plt.show()


##################################################################################
####### Karch algo : CCO####################

np.random.seed(42)

timing = True
writing = False
NTerm = 250

if timing:
    debut = time.time()
    
if writing:
    fd = open('./Results/CCO_output.txt','w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it    
    sys.stdout = fd # Now your file is used by print as destination 

if True:
    #### Parameters to define: ##
    ## About tree
    
    Q_perf = 8.33e3
    N_term = NTerm
    Q_term = Q_perf / N_term
    P_drop = 1.33e7 - 8.38e6 # when Nterm = 4 000, the P_drop is 1.33e7 -7.98e6 #when =Nterm=250 :1.33e7 - 8.38e6
    viscosity = 3.6 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 (check works with radius and length in mm)
    N_con = 20
    N_con_max = 40
    
    # About  convexe perfusion surface : defines a disc surface 
    sphere_center = np.array([50.,50.,50.])
    sphere_radius = 50
    vol_descptr = [sphere_center, sphere_radius]
    volume = np.pi*(sphere_radius**3)*4/3.   
    
    
    #### initialization ##    
    # tree
    #tree_stored= []
    
    #store_cet = []
    tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity)
     
    # source point : define the root position
    root_position = np.array([50,50,100])
    root_node = nclass.Node(0,root_position, Q_perf, -1)
    root_node.set_child_0_index(1)
    root_node.set_label(0)
    tree.add_node(root_node)    
    
    #first segment end: randomly picked inside perfusion surface
    first_node_position = cco_3df.first_segmt_end(vol_descptr)
    
    first_node = nclass.Node(1, first_node_position, Q_perf,0)
    first_node.set_label(0)
    tree.add_node(first_node)
            
                 
    while tree.get_k_term() < N_term: 
        
        success, new_child_location, d_tresh = cco_3df.get_new_location(tree, vol_descptr, N_term)
        
        if (success == False):
            print "impossible to satisfy distance criteria", "d_tresh", d_tresh
            break
        print "new location found"
        cet = []
        adding_location = False
        added_location = []
        
        dtype_r=[("volume", float), ("betasold", float), ("betanew",float), ("branching_locationx",float),("branching_locationy", float) , ("branching_locationz", float),("old_child_index", int)]
        test_N_con_max = False
        
        # test closest neighbors
        neighbors = tree.find_neighbors(new_child_location, N_con)
        for neighbor in neighbors:
            tree_copy = copy.deepcopy(tree)
            convg, res = tree_copy.test_connection(neighbor, new_child_location)    
            if (convg):
                s=res[0]
                cet.append(((s[0]), (s[1][0]), (s[1][1]), (s[2][0]), (s[2][1]),(s[2][2]),(neighbor))) #volume, beta_old, beta_new, branching_location,old_child_index 
        if (len(cet) > 1) or (len(cet)>0 and (tree.get_k_term() == 1)) :
            cet_array = np.array(cet, dtype_r)
            cet_sorted = np.sort(cet_array, order = "volume")
            adding_location = True
            added_location.append(cet_sorted[0])

        # test extra neighbors if no connection candidate has fullfilled constrains
        else: 
            test_N_con_max = False
            neighbors = tree.find_neighbors(new_child_location, N_con_max)
            extra_neighbs = neighbors[N_con:N_con_max]
            for neighbor in extra_neighbs:
                tree_copy = copy.deepcopy(tree)
                convg, res = tree_copy.test_connection(neighbor, new_child_location)    
                if (convg):
                    s=res[0]
                    cet.append(((s[0]), (s[1][0]), (s[1][1]), (s[2][0]), (s[2][1]),(s[2][2]), (neighbor))) #volume, beta_old, beta_new, branching_location,old_child_index
            if (len(cet)) > 1 or (len(cet)>0 and (tree.get_k_term() == 1)):
                print "cet size after Nconmax iter", len(cet)
                cet_array = np.array(cet, dtype_r)
                cet_sorted = np.sort(cet_array, order = "volume")
                adding_location = True
                added_location.append(cet_sorted[0])
        
        if (adding_location): # optimal connection found!
            #store_cet.append(cet)
            if (test_N_con_max):
                print "N con max was tested"
            print "size of cet", len(cet)
            print "CET table", cet
            print "optimal connection from cet ", added_location[-1]
            opt = added_location[-1]
            if (tree.add_connection(opt[6], new_child_location, np.array([opt[3], opt[4], opt[5]]), np.array([opt[1], opt[2]]))):
                branching_added_node = tree.nodes[len(tree.nodes) - 2]       
                added_node = tree.nodes[len(tree.nodes) - 1]
                label = int(len(tree.nodes) / 100.)             
                if (label == 5):
                    label = 4
                added_node.set_label(label)
                branching_added_node.set_label(label)
                print "connection added on tree up node", opt[6]
                print "d_tresh value", d_tresh
                print "k term is now ", tree.get_k_term()
                tree.printing_full()
                #tree_stored.append(copy.deepcopy(tree))
                last_tree = copy.deepcopy(tree)
            else:
                print "failed to add connection on tree"
        else:
            print "location doesn't provide an optimal connection, testing new location"
                
        #keep going until reach Nterm!
        #print "stored cet", store_cet
    fac= 25.
    plot_tree(last_tree, vol_descptr, "./Results/tree_%iNterm%ifac" %(NTerm,fac), fac)#tree_stored[-1]
     




#sys.stdout=old_stdout # here we restore the default behavior
#fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"





    

    
