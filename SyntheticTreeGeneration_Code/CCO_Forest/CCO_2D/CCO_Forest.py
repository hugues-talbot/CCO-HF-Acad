# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import ForestClass as fclass
import CCO_2DFunctions as cco_2df
from pylab import Line2D #figure, gca, 
#import cPickle as pickle
import sys  # Need to have acces to sys.stdout
import time
from multiprocessing import Pool#, TimeoutError
import pickle
import os
import os.path


############# Visualisation tools ####################

def plot_forest(forest, area_descptr, name,w):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    #colors=['k', 'b', 'g', 'm', 'r']
    # labels for cet ['r','m','g','b']
    N=500
    colors = [(1.0,1.0,1.0)]
    colors.extend(plt.cm.jet(np.linspace(0., 1., N)))
    colors.extend([(1.0,1.0,1.0)])
    cmap = mpl.colors.ListedColormap(colors)
    ax.imshow(w,cmap=cmap, vmin=0., vmax=1.)
    ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[1], color = 'r', fill = False))
    ax.add_patch(plt.Circle(area_descptr[0], radius=area_descptr[2], color = 'b', fill = False))
    all_radii = []
    leaves_radii = []
    inv_length_fac = 1. / forest.length_factor
    for tree in forest.trees:
        #print "tree", tree.tree_index
        for sgmt in tree.nodes:
            #print sgmt.index
            #print "parent", sgmt.parent()
            if (sgmt.parent() >= 0):
                #print "segment",sgmt.index
                distal = sgmt.coord
                proximal = tree.get_node(sgmt.parent()).coord
                radius = tree.get_radius(sgmt.index) *inv_length_fac #*10
                all_radii.append(radius)
                if sgmt.is_leaf():
                    leaves_radii.append(radius)               
                ax.add_line(Line2D([distal[0], proximal[0]],[distal[1], proximal[1]], linewidth = radius*2., color = 'k'))#colors[sgmt.label]

    ax.set_xlim(area_descptr[0][0] - area_descptr[1]*1.5,area_descptr[0][0]+ area_descptr[1]*1.5)
    ax.set_ylim(area_descptr[0][1] - area_descptr[1]*1.5,area_descptr[0][1]+ area_descptr[1]*1.5)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    fig.savefig(name+".png")
    fig.savefig(name+".pdf")	
    fig.show()


#######################################################

def test_connection_list(list_input):
    return list_input[0].test_forest_connection(list_input[1], list_input[2],list_input[3])


#########################################################    
    
    
##################################################################################
############# Karch algo : CCO ####################

timing = True
store_data = False
parallelized = True

if timing:
    debut = time.time()
    print debut
if store_data:
    fd = open('./Results/CCO_debug_t.txt','w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it    
    sys.stdout = fd # Now your file is used by print as destination 
    

if True:

    NTerm = 250
    seed = 42

    np.random.seed(seed)
    process_nb = 16
 
    #### Parameters to define: ##
    ## About tree
      
    Q_perf = 8.33e3
    N_term = NTerm
    Q_term = Q_perf / N_term
    gamma = 3.0

    P_drop = 1.33e4 - 8.38e3 # when Nterm = 4 000, the P_drop is 1.33e7 -7.98e6 #when =Nterm=250 :1.33e7 - 8.38e6
    viscosity = 3.6e-3 # 3.6cp = 3.6mPa = 3.6 kg mm-1 s-2 =3.6e-3 Pa.s = 3.6e-9 MPa.s 
    
    N_con = 20
    N_con_max = 40
    
    # About  convexe perfusion surface : defines a disc surface 
    area_center = np.array([120.,130.])#np.array([100.,80.])#
    area_ext_radius =80#100.#50.#
    area_int_radius =37#25#15.#
    area_descptr = [area_center, area_ext_radius, area_int_radius]
    area = np.pi*(area_ext_radius**2 - area_int_radius**2)     
    #potential = cco_2df.potential_image(area_center, area_ext_radius,area_int_radius)


    filename = "potential_rext%i_rint_%i" %(int(area_ext_radius), int(area_int_radius))
    filepath = "./Results/"+filename+".npy"
    if os.path.isfile(filepath):
        print "loading potential from numpy file %s" %filepath
        potential = np.load(filepath)
    else:
        potential = cco_2df.potential_image(area_center, area_ext_radius,area_int_radius)
        np.save("./Results/"+filename, potential)    
    
    
    
    #### initialization ##    
    store_cet = []
    #forest
    forest = fclass.Forest([],NTerm*2,Q_perf*2,P_drop,viscosity,area,np.sqrt(area_ext_radius**2 - area_int_radius**2),potential,area_int_radius,area_center,area_ext_radius,gamma)
    print "forest created"
    # source points : define the root positions
    #first tree
    root_position = np.array([area_center[0],area_center[1]+area_ext_radius])
    first_tree = forest.create_tree(root_position, forest.final_q_perf  * 0.75)#)
    forest.add_tree(first_tree)
    #second tree
    root_position_2 = np.array([area_center[0],area_center[1]-area_ext_radius])
    second_tree = forest.create_tree(root_position_2, forest.final_q_perf  * 0.25)#2. / 3.)
    forest.add_tree(second_tree)
    #third tree
#    root_position_3 = np.array([area_center[0]-area_ext_radius,area_center[1]])
#    third_tree = forest.create_tree(root_position_3,forest.final_q_perf * 2. / 3.)
#    forest.add_tree(third_tree)
    #first segment of each source
    
    print "sources added"
    forest.first_segment_end()
    print "first segmt added"    

    count_extend_neighb_research = np.zeros(3)
    counter = np.zeros(6)
    dead_end_counter = 0
    d_tresh_factor = 1
    while forest.get_fk_term() < N_term*2:

        kterm = forest.get_fk_term()
        success, new_child_location, d_tresh = forest.get_new_location( d_tresh_factor)
        if (success == False):
            print "impossible to satisfy distance criteria", "d_tresh", d_tresh
            break   
        
        cet = [] 
        if parallelized == False:
            cet = [[] for i in range (N_con_max)]
        adding_location = False
        added_location = []
        
        dtype_r=[("convgce", int), ("volume", float), ("betas_and_branching_location", float, (2, 2)),("tree_index", int),("old_child_index", int), ("tree_volume", float)]
                
        test_N_con_max = False
        
        # test closest neighbors
        neighbors = forest.find_forest_neighbors(new_child_location, N_con)
        args = [[forest, neighbors[i][0], neighbors[i][1],new_child_location] for i in range (len(neighbors))]

        #process all neighbors connection test batch by batch
        count_extend_neighb_research[0] = count_extend_neighb_research[0] + 1

        forest.update_forest_length_factor()
        if parallelized == True:
            while (len(cet) < N_con) and (len(cet) < len(neighbors)):           
                end_lim = len(cet) + process_nb if (len(cet) + process_nb < len(neighbors)) else len(neighbors)            
                pool = Pool(processes =  process_nb)               
                res = pool.map(test_connection_list,args[len(cet): end_lim])    
                cet = cet + [res[0]]

                pool.close()
        else:

            for n_index in range(len(neighbors)): 
                print "testing neighbor ", n_index, " ", neighbors[n_index]
                    #break
                res= forest.test_forest_connection(neighbors[n_index][0], neighbors[n_index][1], new_child_location)
                print "res", res  , " kterm ", kterm
                cet[n_index] = res

          
        cet_filtered = filter(None,cet)
        #print "cet_filtered", cet_filtered
        cet_values = np.array(cet_filtered, dtype_r)
        
        #if there is a candidate that converges
        if (np.sum(cet_values['convgce']) >= 1) :        
            cet_sel = cet_values[cet_values['convgce']>0]
            cet_sorted = np.sort(cet_sel, order = "volume")
            cet_final=cet_sorted[0]
            adding_location = True
            added_location.append(cet_final.tolist())
            print "connection to add is", cet_final.tolist()
           
        if (adding_location): # optimal connection found!
            store_cet.append(filter(None,cet))
            opt = added_location[-1]
            #ante_tree = copy.deepcopy(tree)
            
            if (forest.add_connection(opt[2], opt[3], opt[4],new_child_location, opt[5])):
                print "k termmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm is now ", forest.get_fk_term()
                kterm=forest.get_fk_term()
#                if kterm % 10 == 0:
#                    plot_tree(tree, area_descptr, "./Results/InterTree_Nt%i_kt%i_s%i_41d" %(NTerm,kterm,seed),potential)
                #if kterm >600:

                if kterm%50 == 0:
                    plot_forest(forest, area_descptr, "./Results/InterForest_Nt%i_kt%i_s%i_polytree%i_dl75" %(forest.n_term,kterm,seed,len(forest.trees)),potential) 
                if kterm%100 == 0:                    
                    pickle.dump(forest, open("./Results/InterForest_Nt%i_kt_%i_s%i_polytree%i_dl75.p"%(forest.n_term,kterm,seed,len(forest.trees)), "wb"))
#                if kterm == 600:
#                    print "over 600"
#                if kterm == 4:
#                    break  

            else:
                print "failed to add connection on tree"
        else:

            print "ktemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm", forest.get_fk_term()
            print "location doesn't provide an optimal connection, testing new location"
            dead_end_counter = dead_end_counter + 1
            if dead_end_counter == 50:
                dead_end_counter = 0
                d_tresh_factor = 0.8
                print "dead end: decrease d_tresh of 20% to look for new location"

        #keep going until reach Nterm!

    plot_forest(forest, area_descptr, "./Results/Forest_Nt%i_s%i_polytree%i_dl75" %(forest.get_fk_term(),seed,len(forest.trees)), potential)
    pickle.dump(forest, open("./Results/Forest_Nt%i_s%i_polytree%i_dl75.p"%(forest.get_fk_term(),seed,len(forest.trees)), "wb"))



if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"
    

