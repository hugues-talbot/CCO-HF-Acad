# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:35:19 2016

@author: jaquetc
"""

import numpy as np
import Kamiya3D as kami
import NodeClass3D as nclass
import ConvertFunctions as cv
import CCO_3DFunctions as cco_3df
import ForestClass as fclass
#import JsonWriter as jw
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
import copy

############################################################
### LOADING INPUTS

# Location of sources with flow and pressure -> obtained from centerline, mssing branched statistics
sources = np.load("./Inputs/SourcesCor.npy") 
# dtype_r=[("WorldCoordX", float),("WorldCoordY", float), ("WorldCoordZ", float),
#          ("Diameter",float), ("Pressure", float),("Flow", float), 
#          ("VoxelCoordX", float),("VoxelCoordY", float), ("VoxelCoordZ", float)]
#pressure is provided in g mm-1 s-2 = kg m-1 s-2 = Pa
#flow is provided in mm3 s-1

Q_perf = np.sum(sources["Flow"])
viscosity = 3.6e-3 # (Pa s) 
gamma = 3.0    


# matrix to convert between world and voxel coordinates
model_matrix = np.load("./Inputs/ImMatrix.npy")
inv_model_matrix = np.load("./Inputs/ImInvMatrix.npy")
vox_size = np.array([model_matrix[0][0],model_matrix[1][1],model_matrix[2][2]])
vox_volume = vox_size[0]*vox_size[1]*vox_size[2]
#should test if find correspondance between 

# Whole segmented heart-LV potential -> obtained from LargeStructureSegmentationTool
heart_potential = np.load("./Inputs/Heart_potential.npy")
# LV inner vs outer surface potential
lv_potential = np.load("./Inputs/LV_potential.npy")

lv_max_curvature_r = 9. # empirical estimation from mesh measures

vox_perf=np.array(np.where(np.logical_and(lv_potential>0.0, lv_potential < 1.0)))
lv_volume = vox_perf.shape[1] * vox_volume

# Segmented vessels distance map 
segm_vessel_dist = np.load("./Inputs/DistMapFromCenterlines.npy")


##CHECK ALL IMAGE SIZES
im_size_max = np.max(heart_potential.shape)

## CHECK ALL SOURCES ARE INSIDE HEART POT/LV POT ==>otherwise need to change potential
## if outside LV: calculate dist max and see if sources are too far from LV --> use projection
##################################################################################
############# INITIALIZATION ####################

timing = True
store_data = True
parallelized = False
filename = "./Results/Test"

if timing:
    debut = time.time()
    print debut
if store_data:
    fd = open(filename + ".txt",'w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it     
    sys.stdout = fd # Now your file is used by print as destination 
    

if True:
  
    seed = 42
    np.random.seed(seed)
    process_nb = 16
 
    #### Parameters to define: ##
    ## About tree
    NTerm = 250 
    InterTerm = 0
    P_term = 8.38e3 #(Pa)
    Q_term = Q_perf / NTerm
    N_con = 20

    r_f = np.power(3*lv_volume /(4*np.pi),1./3)
    forest = fclass.Forest([], NTerm, Q_perf, P_term, viscosity, lv_volume, r_f, vox_size, heart_potential, lv_potential,segm_vessel_dist, lv_max_curvature_r, lv_potential.shape, gamma)
    print "forest created"
    #check all sources location are inside perf territory
    #then add them to forest
    for source in sources:
        forest.create_tree(np.array([source["VoxelCoordZ"], source["VoxelCoordY"], source["VoxelCoordX"]]), source["Flow"], source["Pressure"])

    print "trees added"
    #first segment end: randomly picked inside perfusion surface
    if InterTerm >0:          
        surface = True
        surface_tol = 0.75 #from inside lv?
    else: 
        surface = False
        surface_tol = 0.
            
    
    forest.first_segment_end(surface, surface_tol) #need to update the ffunction when will do stage growth
    print "first segment assigned"
    dead_end_counter = 0
    d_tresh_factor = 1
         
    while forest.get_fk_term() < NTerm:
        kterm = forest.get_fk_term()
            
        if kterm < InterTerm:          
            surface = True
        else: 
            surface = False
            surface_tol = 0.
        
        success, new_child_location, d_tresh = forest.get_new_location(d_tresh_factor)
        if (success == False):
            print "impossible to satisfy distance criteria", "d_tresh", d_tresh
            break       
        
        cet = [] 
        if parallelized == False:
            cet = [[] for i in range (N_con)]
        adding_location = False
        added_location = []
        
        dtype_r=[("convgce", int), ("volume", float), ("betas", float,(2)),("branching_location",float,(3)),("tree_index", int),("old_child_index", int), ("tree_volume", float)]
                
        test_N_con_max = False
        
        # test closest neighbors
        neighbors = forest.find_forest_neighbors(new_child_location, N_con)


        args = [[forest, neighbors[i][0], neighbors[i][1],new_child_location, surface, surface_tol] for i in range (len(neighbors))]

        #process all neighbors connection test batch by batch
        if parallelized == True:
            while (len(cet) < N_con) and (len(cet) < len(neighbors)):           
                end_lim = len(cet) + process_nb if (len(cet) + process_nb < len(neighbors)) else len(neighbors)            
                pool = Pool(processes =  process_nb)               
                res = pool.map(cco_3df.test_connection_list,args[len(cet): end_lim])    
                cet = cet + [res[0][1:]]
                pool.close()
        else:
            for n_index in range(len(neighbors)):
                res= forest.test_forest_connection(neighbors[n_index][0], neighbors[n_index][1], new_child_location, surface, surface_tol)                 
                cet[n_index] = res
                print "res", res                 

        
        cet_filtered = filter(None,cet)
        cet_values = np.array(cet_filtered, dtype_r)
        
        #if there is a candidate that converges
        if (np.sum(cet_values['convgce']) >= 1) or (np.sum(cet_values['convgce']) > 0 and forest.get_fk_term() == 1):        
            cet_sel = cet_values[cet_values['convgce']>0]
            cet_sorted = np.sort(cet_sel, order = "volume")
            cet_final=cet_sorted[0]
            adding_location = True
            added_location.append(cet_final.tolist())
           
        if (adding_location): # optimal connection found!
            opt = added_location[-1]   
            if (forest.add_connection([opt[2], opt[3]], opt[4], opt[5],new_child_location, opt[6])):
                print "k termmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm is now ", forest.get_fk_term()
                kterm=forest.get_fk_term()
                d_tresh_factor = 1.
                if kterm == 50:
                    break
#                if kterm%10 == 0:
#                    name =filename+"_F_Nt%i_kt%i_s%i_ellip" %(NTerm,kterm,seed)
#                    pickle.dump(forest, open(name + ".p", "wb"))
            else:
                print "failed to add connection on tree"
        else:              
            print "ktemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm", forest.get_fk_term()
            print "location doesn't provide an optimal connection, testing new location"
            dead_end_counter = dead_end_counter +1
            print "dead_end_counter = ", dead_end_counter
            if dead_end_counter == 10: #if too small doesn't produce homogeneous tree?
                dead_end_counter = 0
                d_tresh_factor = d_tresh_factor *0.85
                print "dead end: decrease d_tresh of 20% to look for new location", d_tresh_factor

        #keep going until reach Nterm!
    
    
    print "CCO done"
    name = filename+"_Nt%i_kt%i_s%i" %(NTerm, forest.get_fk_term(),seed)  
    #cco_3df.plot_forest(forest, name)
    #pickle.dump(forest, open(name + ".p", "wb"))
    #pickle.dump(forest.trees,open(name + ".p", "wb"))
    pickle.dump([copy.copy(i) for i in forest.trees],open(name + ".p", "wb"))
    #print "forest saved"
    cv.write_json(forest,model_matrix,"C:/Users/cjaquet/Documents/SynthTreeData/3c9e679d-2eab-480c-acaa-31da12301b0a/Forest.json")
    #print "json generated"
    
if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin =time.time()
    print "duration = ", fin-debut, "secondes"
    

