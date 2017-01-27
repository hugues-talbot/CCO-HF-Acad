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
#from pylab import figure, gca, Line2D
#import pylab as pl
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
#for each input should make test: if npy file exist load it otherwise create data and save them

# Location of sources with flow and pressure -> obtained from centerline, mssing branched statistics
print "loading sources location"
source_path = "./Inputs/Sources_LADLCX.npy"
if os.path.isfile(source_path):
    sources = np.load(source_path)
else:
    print "generate sources"
    sources = cv.get_data_fp("./Inputs/SourcesCorrectedWithDiam.txt") #AllSourcesInit, #LADLCXPressureAndFlowInputs
    np.save(source_path, sources)
       
# dtype_r=[("WorldCoordX", float),("WorldCoordY", float), ("WorldCoordZ", float),
#          ("Diameter",float), ("Pressure", float),("Flow", float), 
#          ("VoxelCoordX", float),("VoxelCoordY", float), ("VoxelCoordZ", float)]
#pressure is provided in g mm-1 s-2 = kg m-1 s-2 = Pa
#flow is provided in mm3 s-1

Q_perf = np.sum(sources["Flow"])
viscosity = 3.6e-3 # (Pa s) 
gamma = 3.0    


# matrix to convert between world and voxel coordinates
print "loading matrices"
model_matrix = np.load("./Inputs/ImMatrix.npy")
inv_model_matrix = np.load("./Inputs/ImInvMatrix.npy")
vox_size = np.array([model_matrix[0][0],model_matrix[1][1],model_matrix[2][2]])
vox_volume = vox_size[0]*vox_size[1]*vox_size[2]
#should test if find correspondance between 

# Whole segmented heart-LV potential -> obtained from LargeStructureSegmentationTool
print "loading potential images"

hp_path = "./Inputs/heart_potential.npy"
if os.path.isfile(hp_path):
    heart_potential = np.load(hp_path)
else:
    print "generate heart potential"
    heart_potential = cv.generate_heart_potential("./Inputs/heart.mha","./Inputs/lv.mha")
    np.save(hp_path, heart_potential)
    
la_path = "./Inputs/l_a.npy"
if os.path.isfile(la_path):
    left_a = np.load(la_path)
else:
    print "load left atrium mha"
    left_a = cv.open_mha_mha("./Inputs/left_atrium.mha")
    np.save(la_path, left_a)
#heart_potential[np.where(left_a==1)]=2. #to avoid segments crossing at top of concavity
    
# LV inner vs outer surface potential
lv_path = "./Inputs/LV_potential.npy"
if os.path.isfile(lv_path):
    lv_potential = np.load(lv_path)
else:
    print "generate lv potential"
    lv_potential = cv.generate_lv_potential("./Inputs/lv.mha")
    np.save(lv_path, lv_potential)
lv_pot_zero = np.zeros(lv_potential.shape) + lv_potential
lv_pot_zero[np.where(lv_pot_zero<0.)] = 0.
lv_pot_zero[np.where(lv_pot_zero>0.999999)] = 2.
lv_max_curvature_r = 9. # empirical estimation from mesh measures

vox_perf=np.array(np.where(np.logical_and(lv_potential>0.0, lv_potential < 1.0)))
lv_volume = vox_perf.shape[1] * vox_volume

# Segmented vessels distance map 
print "loading distance map"
sgm_vess_path = "./Inputs/DistMapFromCenterlines.npy"
if os.path.isfile(sgm_vess_path):
    segm_vessel_dist = np.load(sgm_vess_path)
else:
    segm_vessel_dist = cv.open_mha("./Inputs/fmm_path_with_LCX_and_proj_and_inner5_and3.mha")
    np.save(sgm_vess_path, segm_vessel_dist)

##CHECK ALL IMAGE SIZES
im_size_max = np.max(heart_potential.shape)

## CHECK ALL SOURCES ARE INSIDE HEART POT/LV POT ==>otherwise need to change potential
## if outside LV: calculate dist max and see if sources are too far from LV --> use projection
##################################################################################
############# INITIALIZATION ####################

timing = True
store_data = True
parallelized = False
generate_obj = True
filename = "./Results/Dbg_f"
name_filemaya = "Dbg_f"
path_out = "C:/Users/cjaquet/Documents/SynthTreeData/3c9e679d-2eab-480c-acaa-31da12301b0a/ResultsForMaya/"
ktermbreak = len(sources) + 196
NTerm = 500 
InterTerm = len(sources) + 196

if timing:
    debut = time.time()
    print debut
if store_data:
    fd = open(filename + ".txt",'w') # open the result file in write mode
    old_stdout = sys.stdout   # store the default system handler to be able to restore it     
    sys.stdout = fd # Now your file is used by print as destination 
    

if False:
    #def regenerate_forest

    P_term = 8.38e3 #(Pa) 
    Q_term = Q_perf / NTerm
    r_f = np.power(3*lv_volume /(4*np.pi),1./3)
    forest = fclass.Forest(pickle.load(open("./Results/NewVersions_Nt500_kt79_s42.p")), NTerm, Q_perf, P_term, viscosity, lv_volume, r_f, vox_size, heart_potential, left_a, lv_pot_zero,segm_vessel_dist, lv_max_curvature_r, lv_potential.shape, gamma)
    forest.update_forest_length_factor()    
    cv.write_json(forest,model_matrix,path_out +"ForestNewRadius.json")



if True:
    
    seed = 42
    np.random.seed(seed)
    process_nb = 16
 
    #### Parameters to define: ##
    ## About tree

    P_term = 8.38e3 #(Pa)
    Q_term = Q_perf / NTerm
    N_con = 20

    r_f = np.power(3*lv_volume /(4*np.pi),1./3)
    forest = fclass.Forest([], NTerm, Q_perf, P_term, viscosity, lv_volume, r_f, vox_size, heart_potential, left_a, lv_pot_zero,segm_vessel_dist, lv_max_curvature_r, lv_potential.shape, gamma)
    print "forest created"
    #check all sources location are inside perf territory
    # and at a max distance to lv:
    
        
    
    
    #then add them to forest
    ind = 0
    source_to_be_deleted=[]
    for source in sources: #check if they are inside heart or lv
        source_coord=np.array([source["VoxelCoordZ"], source["VoxelCoordY"], source["VoxelCoordX"]])
        source_vect_coord=  np.array([source["VectVoxelCoordZ"], source["VectVoxelCoordY"], source["VectVoxelCoordX"]])
        ins=True
#        if ind >0 :
#             ins, val = forest.inside_heart(np.array([source["VoxelCoordZ"], source["VoxelCoordY"], source["VoxelCoordX"]]))
#        if (ins==True or (val > 0.999)):
#             forest.create_tree(np.array([source["VoxelCoordZ"], source["VoxelCoordY"], source["VoxelCoordX"]]), np.array([source["VectVoxelCoordZ"], source["VectVoxelCoordY"], source["VectVoxelCoordX"]]), source["Flow"], source["Pressure"], source["Diameter"])
#             ind = ind + 1 
#        else:
#            print "source out of heart and lv, need to update the potential to integrate it", val
#            break
        ins, val = forest.inside_heart(source_coord)
#        if ind > ktermbreak:
#            print "out source", ind
#            break
        if (val >=1):
            forest.create_tree(source_coord,source_vect_coord, source["Flow"], source["Pressure"], source["Diameter"])
            ind = ind + 1 
        else:
            if ins==True :
                #test dist to lv
                end= forest.short_segmt_end(source_coord, 50., True,np.zeros(3))
                closest_to_lv = forest.dist_to_lv_via_sampling(source_coord, end, 100)
                dist_lv = cco_3df.length(closest_to_lv - source_coord, vox_size)
                print "ind", ind,"dist_lv", dist_lv
                if dist_lv > 15:
                    source_to_be_deleted.append(ind)
                    ind = ind+1
                else:
                    forest.create_tree(source_coord,source_vect_coord, source["Flow"], source["Pressure"], source["Diameter"])
                    ind = ind+1
            else:
                print "source out of heart and lv, need to update the potential to integrate it", val
                source_to_be_deleted.append(ind)
                ind = ind+1
                break
        
    print "trees added", len(forest.trees)
    print "source to be deleted", source_to_be_deleted
    #first segment end: randomly picked inside perfusion surface
    if InterTerm >0:          
        surface = True
        surface_tol = 1.5 #value in (heart+lv) potential
    else: 
        surface = False
        surface_tol = -1.
            
    
    forest.first_segment_end(surface_tol) #need to update the ffunction when will do stage growth
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
        
        success, new_child_location, d_tresh = forest.get_new_location(d_tresh_factor, surface_tol)
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
                print "n_index", n_index, len(neighbors), N_con                
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
            print "cet sorted", cet_sorted 
            opt = added_location[-1]   
            if (forest.add_connection([opt[2], opt[3]], opt[4], opt[5],new_child_location, opt[6])):
                print "k termmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm is now ", forest.get_fk_term()
                kterm=forest.get_fk_term()
                d_tresh_factor = 1.
                if (kterm == ktermbreak):
                    break
                if kterm%10 == 0 :
#                if kterm > 70 :
                    cv.write_json(forest,model_matrix,path_out + "DbgF%i.json" %kterm)
                if kterm == InterTerm:
                    cv.write_json(forest,model_matrix,path_out + "DbgF%i.json" %kterm)                    
#                if kterm>70 and kterm%5==0:
#                    cv.write_json(forest,model_matrix,path_out + "SpForest%i.json" %kterm)
##                    name =filename+"_F_Nt%i_kt%i_s%i_ellip" %(NTerm,kterm,seed)
#                    pickle.dump(forest, open(name + ".p", "wb"))
            else:
                print "failed to add connection on tree"
        else:              
            print "ktemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm", forest.get_fk_term()
            print "location doesn't provide an optimal connection, testing new location"
            dead_end_counter = dead_end_counter +1
            print "dead_end_counter = ", dead_end_counter
            if dead_end_counter == 5: #if too small doesn't produce homogeneous tree?
                dead_end_counter = 0
                d_tresh_factor = d_tresh_factor *0.85
                print "dead end: decrease d_tresh of 20% to look for new location", d_tresh_factor

        #keep going until reach Nterm!
#    for tree in forest.trees:
#        c = tree.nodes[0].coord
#        d= tree.nodes[1].coord
#        print tree.tree_index, cco_3df.length(c-d, vox_size), forest.dist_max[tree.tree_index]
#    for tree in forest.trees:
#        print tree.tree_index,  tree.get_h(tree.nodes[0].coord),tree.get_h(tree.nodes[1].coord), "length",tree.length(1)
    print "CCO done"
    forest.printing()
    name = filename+"_Nt%i_kt%i_s%i" %(NTerm, forest.get_fk_term(),seed)  
    #cco_3df.plot_forest(forest, name)
    #pickle.dump(forest, open(name + ".p", "wb"))
    #pickle.dump(forest.trees,open(name + ".p", "wb"))
    pickle.dump([copy.copy(i) for i in forest.trees],open(name + ".p", "wb"))
    #print "forest saved"
    cv.write_json(forest,model_matrix,path_out +name_filemaya+".json")
    #print "json generated"
    #print "sources", sources
if store_data:
    sys.stdout=old_stdout # here we restore the default behavior
    fd.close() # to not forget to close your file

if timing:
    fin = time.time()
    print "duration = ", fin-debut, "secondes"
    

if generate_obj:
    os.chdir("C:\\home\\GitMod\\hfm\\install_build_Release\\bin")
    seg_tool = "C:\\home\\GitMod\\hfm\\build_Release\\bin\\Release\\SyntheticTreeGeneration.exe"
    input_zhf = "C:\\Users\\cjaquet\\Documents\\SynthTreeData\\3c9e679d-2eab-480c-acaa-31da12301b0a\\ResultsForMaya\\"+  name_filemaya+".json"
    output_dir = "C:\\Users\\cjaquet\\Documents\\SynthTreeData\\3c9e679d-2eab-480c-acaa-31da12301b0a"
    seg_args = "--input_zhf {0} --output_directory {1} ".format(input_zhf, output_dir)
    print seg_args
    os.system(seg_tool + " " + seg_args)
    os.chdir("C:\\Users\\cjaquet\\Documents\\GitAcad\\HeartFlow\\SyntheticTreeGeneration_Code\\Patient_CCO")