# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import Kamiya3D as kami
import NodeClass3D as nclass
import CCO_3DFunctions as cco_3df
import copy
from pylab import figure, gca, Line2D
import pylab as pl
import random
import os
import os.path
import cPickle as pickle
import sys  # Need to have acces to sys.stdout
import time
from multiprocessing import Pool, TimeoutError
import pickle
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d



############# Visualisation tools ####################

def plot_tree(tree, vol_descptr, name, factor_radius):
    wid=16 # Must be proportional to x and y limits below
    hei=16
    ax = a3.Axes3D(pl.figure(figsize=(wid, hei)))
    
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = vol_descptr[0][0]+ vol_descptr[1] * np.outer(np.cos(u), np.sin(v))
    y = vol_descptr[0][1]+ vol_descptr[1] * np.outer(np.sin(u), np.sin(v))
    z = vol_descptr[0][2]+ vol_descptr[1] * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1) #,
    
    x = vol_descptr[0][0]+ vol_descptr[2] * np.outer(np.cos(u), np.sin(v))
    y = vol_descptr[0][1]+ vol_descptr[2] * np.outer(np.sin(u), np.sin(v))
    z = vol_descptr[0][2]+ vol_descptr[2] * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z,rstride=2, cstride=2, linewidth=0,color='r',alpha = 0.1)
    
    #setting figure so that we get linewidth in data unit
    bound = vol_descptr[1]*1.5
    ax.set_xlim(vol_descptr[0][0] - bound,vol_descptr[0][0]+ bound)#using same bound diff for all axis 
    ax.set_ylim(vol_descptr[0][1] - bound,vol_descptr[0][1]+ bound)#--> proportional to the figure size
    ax.set_zlim(vol_descptr[0][2] - bound,vol_descptr[0][2]+ bound)
    ax.axes.set_aspect('equal')
    x1,x2,y1,y2 = ax.axis()
    print "axis",x1,x2,y1,y2
    yrange =   y2 - y1 #don't know if necessary to calculate pointlinewid_factor
    # line is in points: 72 points per inch
    point_hei=hei*72 
    # For the calculation below, you have to adjust width by 0.8
    # because the top and bottom 10% of the figure are labels & axis
    pointlinewid_factor = point_hei * 0.8  # corresponding width in pts ( /yrange ?)
    
    inv_length_fac = 1. / tree.length_factor
    for sgmt in tree.nodes:
        if (sgmt.parent() >= 0):
            distal = sgmt.coord
            proximal = tree.get_node(sgmt.parent()).coord
            radius = tree.get_radius(sgmt.index) *inv_length_fac
            print "radius", radius
            verts = [zip([distal[0], proximal[0]],[distal[1], proximal[1]],[distal[2], proximal[2]])]
            tri = a3.art3d.Poly3DCollection(verts)
            tri.set_color('k')
           

            tri.set_linewidth(radius*2.*pointlinewid_factor/(vol_descptr[1]*3.) )
            ax.add_collection3d(tri)
            


    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.grid(False)
    plt.savefig(name+".png")
    plt.savefig(name+".pdf")
    plt.show()


    
##################################################################################
############# Karch algo : CCO ####################

timing = True
store_data = False
parallelized = False

if timing:
    debut = time.time()
    print debut
if store_data:
    fd = open('./Results/CCO3D_debug.txt','w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it    
    sys.stdout = fd # Now your file is used by print as destination 
    

if True:

    NTerm = 4000
    seed = 42

    np.random.seed(seed)
    process_nb = 16
 
    #### Parameters to define: ##
    ## About tree
      
    Q_perf = 8.33e3
    N_term = NTerm
    Q_term = Q_perf / N_term

    P_drop = 1.33e4 - 7.98e3 # when Nterm = 4 000, the P_drop is 1.33e7 -7.98e6 #when =Nterm=250 :1.33e7 - 8.38e6
    viscosity = 3.6e-3 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 =3.6e-3 Pa.s = 3.6e-9 MPa.s 
    gamma = 3.0    
    # need to be carefull about unit: relativ between P_drop and viscosity
    # ex: 1.33e7  - 7.98e6 (kg mm-1 s-2)   and 3.6    (kg mm-1 s-2)
    # ex: 1.33e-2 - 7.98e-3 (MPa = N mm-2) and 3.6e-9 (MPa)
    # ex: 1.33e4  - 7.98e3 (Pa)            and 3.6e-3 (Pa.s)  
    N_con = 20

    # About  convexe perfusion surface : defines a disc surface 
    v_center = np.array([14.,14., 14.])#np.array([80.,80.,80.])#
    v_ext_radius =10.#50#
    v_int_radius =4.#15#
    v_descptr = [v_center, v_ext_radius, v_int_radius]
    vperf = (4./3.) * np.pi*(v_ext_radius**3 - v_int_radius**3) 
    filename = "potential_rext%i_rint_%i" %(int(v_ext_radius), int(v_int_radius))
    filepath = "./Results/"+filename+".npy"
    if os.path.isfile(filepath):
        potential = np.load(filepath)
    else:
        potential = cco_3df.potential_image(v_center, v_ext_radius,v_int_radius)
        np.save("./Results/"+filename, potential)
    
    #### initialization ##    
    store_cet = []
    tree = nclass.Tree([], N_term, Q_perf, P_drop, viscosity, vperf, np.power((v_ext_radius**3 - v_int_radius**3),1/3.), potential, v_int_radius, v_center, v_ext_radius, gamma)
     
    # source point : define the root position
    root_position = np.array([v_center[0],v_center[1]+v_ext_radius, v_center[2]])
    root_node = nclass.Node(0,root_position, Q_perf, -1)
    root_node.set_child_0_index(1)
    root_node.set_label(0)
    tree.add_node(root_node)    
    #first segment end: randomly picked inside perfusion surface
    tree.update_length_factor()
    first_node_position = tree.first_segmt_end()
    first_node = nclass.Node(1, first_node_position, Q_perf,0)
    first_node.set_label(0)
    tree.add_node(first_node)
    tree.update_flow()
    tree.update_length_factor()
    tree.depthfirst_resistances(0)        

    counter = np.zeros(6)
    dead_end_counter = 0
    d_tresh_factor = 1
    while tree.get_k_term() < N_term: 
        success, new_child_location, d_tresh = tree.get_new_location(N_term, d_tresh_factor)
        if (success == False):
            print "impossible to satisfy distance criteria", "d_tresh", d_tresh
            break       
        
        cet = [] 
        if parallelized == False:
            cet = [[] for i in range (N_con)]
        adding_location = False
        added_location = []
        
        dtype_r=[("convgce", int), ("volume", float), ("betas", float,(2)),("branching_location",float,(3)),("old_child_index", int)]
                
        test_N_con_max = False
        
        # test closest neighbors
        neighbors = tree.find_neighbors(new_child_location, N_con)
        args = [[tree, neighbors[i],new_child_location] for i in range (len(neighbors))]

        #process all neighbors connection test batch by batch
#        count_extend_neighb_research[0] = count_extend_neighb_research[0] + 1
        if parallelized == True:
            while (len(cet) < N_con) and (len(cet) < len(neighbors)):           
                end_lim = len(cet) + process_nb if (len(cet) + process_nb < len(neighbors)) else len(neighbors)            
                pool = Pool(processes =  process_nb)               
                res = pool.map(cco_3df.test_connection_list,args[len(cet): end_lim])    
                cet = cet + [res[0][1:]]
                #getting stat
                code= res[0][0]
                counter[0]= counter[0]+1
                if code>0:
                    if code>2:
                        counter[1] = counter[1]+1
                        counter[2] = counter[2]+ code
                    if code==2:
                        counter[3] = counter[3]+1
                    if code==1:
                        counter[3] = counter[3]+1
                        counter[4] = counter[4]+1
                else:
                    counter[5]=counter[5]+1
                pool.close()
        else:
            for n_index in range(len(neighbors)):
                tree_copy = copy.deepcopy(tree)
                res= tree_copy.test_connection(neighbors[n_index], new_child_location)
                cet[n_index] = res[1:]
                #getting stat
                counter[0]= counter[0]+1
                if res[0]>0:
                    if res[0]>2:
                        counter[1] = counter[1]+1
                        counter[2] = counter[2]+ res[0]
                    if res[0]==2:
                        counter[3] = counter[3]+1
                    if res[0]==1:
                        counter[3] = counter[3]+1
                        counter[4] = counter[4]+1
                else:
                    counter[5]=counter[5]+1
                

        
        cet_filtered = filter(None,cet)
        cet_values = np.array(cet_filtered, dtype_r)
        
        #if there is a candidate that converges
        if (np.sum(cet_values['convgce']) >= 1) or (np.sum(cet_values['convgce']) > 0 and tree.get_k_term() == 1):        
            cet_sel = cet_values[cet_values['convgce']>0]
            cet_sorted = np.sort(cet_sel, order = "volume")
            cet_final=cet_sorted[0]
            adding_location = True
            added_location.append(cet_final.tolist()[1:])
           
        if (adding_location): # optimal connection found!
            store_cet.append(filter(None,cet))
            opt = added_location[-1]
           
            if (tree.add_connection(opt[3], new_child_location, opt[2], opt[1])):
                print "k termmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm is now ", tree.get_k_term()
                kterm=tree.get_k_term()
                if kterm%100 == 0:
                    plot_tree(tree, v_descptr, "./Results/InterTree_Nt%i_kt%i_s%i" %(NTerm,kterm,seed),potential) 

                if kterm == 1500:
                    break          
            else:
                print "failed to add connection on tree"
        else:

            print "ktemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm", tree.get_k_term()
            print "location doesn't provide an optimal connection, testing new location"
            dead_end_counter = dead_end_counter +1
            if dead_end_counter == 50:
                dead_end_counter = 0
                d_tresh_factor = 0.8
                print "dead end: decrease d_tresh of 20% to look for new location"

        #keep going until reach Nterm!
        if tree.get_k_term() == 10:
            break
    plot_tree(tree, v_descptr, "./Results/tree_Nt%i_s%i" %(tree.get_k_term(),seed), potential)
    pickle.dump(tree, open("./Results/treetNt%i_s%i.p"%(tree.get_k_term(),seed), "wb"))



if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"
    

