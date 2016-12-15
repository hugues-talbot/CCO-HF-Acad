# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:32:18 2016

@author: cjaquet
"""

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import CCO_3DFunctions as cc3df

def open_mha(mha_file):    
    img_reader = vtk.vtkMetaImageReader()
    img_reader.SetFileName(mha_file)
    img_reader.Update()
    vtk_marker = img_reader.GetOutput()
    mk_dim = np.array(vtk_marker.GetDimensions())
    mk_rdim = np.array((mk_dim[2], mk_dim[1], mk_dim[0]))
    print "dim", mk_rdim
    np_marker = vtk_to_numpy(vtk_marker.GetPointData().GetScalars()).reshape(mk_rdim)
    #np_marker[np.where(np_marker > 1000)] = 0
    np_marker = np_marker.astype(float)
    return np_marker
    

def voxel_to_world(voxel, matrix):
    w = float(voxel[0]*matrix[0][3] + voxel[1]*matrix[1][3] + voxel[2]*matrix[2][3] + matrix[3][3]) 
    res= np.zeros(len(voxel))
    for i in range(len(res)):
        res[i] = (voxel[0]*matrix[0][i] + voxel[1]*matrix[1][i] + voxel[2]*matrix[2][i] + matrix[3][i]) / w
    return res;
    
def world_to_voxel(coord, inv_matrix):
    res = voxel_to_world(coord, inv_matrix)
    int_res = (res + 0.5).astype(int)    
    return res, int_res
    
def world_to_sdf(world_coord,transl_offset, sdf_voxel_size):
    return (world_coord - transl_offset)/float(sdf_voxel_size)
    
def sdf_to_world(sdf_coord, transl_offset, sdf_voxel_size):
    return sdf_coord*sdf_voxel_size + transl_offset

#cv.get_data_fp("./Inputs/LocationDiameterPressureFlowInputs_LAD.txt")    
def get_data_fp(input_file):
    with open(input_file, "r") as filep:
        content = [x.split('\t') for x in filep.read().splitlines()]
        tb = []
        for j in content: 
            a= [float(i) for i in j]
            tb.append(tuple(a))
    dtype_r=[("WorldCoordX", float),("WorldCoordY", float), ("WorldCoordZ", float),("Diameter",float), ("Pressure", float),("Flow", float), ("VoxelCoordX", float),("VoxelCoordY", float), ("VoxelCoordZ", float),("VectVoxelCoordX", float),("VectVoxelCoordY", float), ("VectVoxelCoordZ", float)]
    array = np.array(tb, dtype = dtype_r )
    array_sorted = np.sort(array, order = "Flow")[::-1]  
    return array_sorted



import json

def write_json(forest, matrix, filename):
    file_out = open(filename, 'w')
    fdict={}
    fdict.update({"Flow" : forest.final_q_perf})
    fdict.update({"Trees" : len(forest.trees)})
    fdict.update({"Nterm" : forest.n_term})
    fdict.update({"FKterm" : forest.get_fk_term()})
    fdict["Trees"]= {}
    for tree in forest.trees:
        fdict["Trees"][tree.tree_index]={}
        fdict["Trees"][tree.tree_index].update({"Tree flow" : tree.final_q_perf})
        fdict["Trees"][tree.tree_index].update({"Tree pressure drop" : tree.p_drop})
        fdict["Trees"][tree.tree_index].update({"Kterm" : tree.get_k_term()})
        fdict["Trees"][tree.tree_index].update({"Resistance": tree.resistance(tree.get_root_index())})
        fdict["Trees"][tree.tree_index]["Nodes"]= {}
        for node in tree.nodes:
            fdict["Trees"][tree.tree_index]["Nodes"][node.index]={}
            print "index", node.index, "node.coord", node.coord    
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Location" : (voxel_to_world(node.coord[::-1], matrix)).tolist()})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Parent" : node.parent_index})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Children" : node.children_index.tolist()})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Radius" : tree.get_radius(node.index)})
            
    json.dump(fdict, file_out, indent = 4, sort_keys = True)
    file_out.close()
    return True




if False:
    N=50
    colors = [(1.0,1.0,1.0)]
    colors.extend(plt.cm.jet(np.linspace(0., 1., N)))
    colors.extend([(1.0,1.0,1.0)])
    cmap =plt.colors.ListedColormap(colors) #plt.cm.jet


if False:
    #load txt file and store all in a structured array
    with open("./Inputs/LocationDiameterPressureFlowInputs_LAD.txt", "r") as filep:
        content = [x.split('\t') for x in filep.read().splitlines()]
        tb = []
        for j in content: 
            a= [float(i) for i in j]
            tb.append(tuple(a))
    dtype_r=[("WorldCoordX", float),("WorldCoordY", float), ("WorldCoordZ", float),("Diameter",float), ("Pressure", float),("Flow", float), ("VoxelCoordX", float),("VoxelCoordY", float), ("VoxelCoordZ", float)]
    
    array = np.array(tb, dtype = dtype_r )
    array_sorted = np.sort(array, order = "Flow")[::-1]  
        
    

import SimpleITK as sitk

def generate_lv_potential(filepath):
    lv = open_mha(filepath)
    ###generate outer marker (smooth lv with light closing)
    lv_int = lv.astype(int)
    lv_im = sitk.GetImageFromArray(lv_int)
    lv_smooth = sitk.BinaryMorphologicalClosing(lv_im, 5)
    lv_smth = sitk.GetArrayFromImage(lv_smooth)

    ###generate inner marker:  light opening(outermarker -heavy closing(lv))
    closing_fac = 100
    rsult = sitk.BinaryMorphologicalClosing(lv_im, closing_fac)
    sub = rsult - lv
    sub_im = sitk.GetImageFromArray(sub.astype(int))
    sub_open = sitk.BinaryMorphologicalOpening(sub_im, 5)
    sub_a = sitk.GetArrayFromImage(sub_open)
    
    #generate lv potential
    out_inv = lv_smth +1
    out_inv[np.where(out_inv>1)]=0
    out_marker = lv_smth[40:200,100:,150:]
    in_marker = sub_a[40:200,100:,150:]
    
    mask = (in_marker + out_marker) +1
    mask[np.where(mask>1)]=0 
    cc3df.generate_potential(mask, inner[40:200,100:,150:], out_inv[40:200,100:,150:])
    dimb = inner.shape
    final_size = np.zeros(dimb)-1
    final_size[40:200,100:,150:] = pot
    return final_size
    
def generate_heart_potential(filepath1, filepath2):
    from skimage.segmentation import random_walker
    from skimage.morphology import binary_dilation, binary_erosion
    lv = open_mha(filepath2)
    heart = open_mha(filepath1)  
    dim_h = lv.shape

    heart_dil = binary_dilation(heart[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50], np.ones((15,15,15)))    
    heart_mask=heart_dil+1
    heart_mask[np.where(heart_mask>1)]=0     
    lv_mask = binary_erosion(lv[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50], np.ones((5,5,5)))
    mask = (lv_mask + heart_mask) 

    
    inner_marker = lv[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50]
    heart_small_dil =  binary_dilation(heart, np.ones((10,10,10)))
    heart_marker = heart_small_dil +1
    heart_marker[np.where(heart_marker>1)]=0 
    outer_marker = heart_marker[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50]
    markers = outer_marker*2 + inner_marker
    
    result = random_walker(mask, markers, copy =True, return_full_prob = True, mode= 'cg_mg')
    
    result[0][np.where(result[0] == 0.)] = -1 
    final_size = np.zeros(dim_h)-1
    final_size[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50] = result[0]
    return final_size
 


   
if False:
    #generate markers for lv potential
    
    
    ####load lv
    lv = open_mha("C:/Users/cjaquet/Documents/SynthTreeData/3c9e679d-2eab-480c-acaa-31da12301b0a/lv.mha")  
    
    ####generate outer marker (smooth lv with light closing)
    #lv_int = lv.astype(int)
    #lv_im = sitk.GetImageFromArray(lv_int)
    #lv_smooth = sitk.BinaryMorphologicalClosing(lv_im, 5)
    #lv_smth = sitk.GetArrayFromImage(lv_smooth)
    #np.save("lv_smooth.npy", lv_smth)
    lv_smth = np.load("lv_outter_marker.npy")
    
    
    ####generate inner marker:  light opening(outermarker -heavy closing(lv))
    #closing_fac = 100
    #r = sitk.BinaryMorphologicalClosing(lv_im, closing_fac)
    rsult = np.load("Closing100.npy")
    #sub = rsult - lv
    #sub_im = sitk.GetImageFromArray(sub.astype(int))
    #sub_open = sitk.BinaryMorphologicalOpening(sub_im, 5)
    #sub_a = sitk.GetArrayFromImage(sub_open)
    #np.save("lv_inner_marker.npy", sub_a)
    sub_a = np.load("lv_inner_marker.npy")
    
    #### sum them for visualization
    lv_reverse = lv_smth+1
    lv_reverse[np.where(lv_reverse>1)]=0
    markers = lv_smth*2 + sub_a
    for i in range (6):    
        a= markers[40 + i*30]
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        ax.imshow(a)
        fig.show()
        fig.savefig("Markers_"+str(40 + i*30)+".png")
      
if False :
  #generate lv potential
  outter = np.load("lv_outter_marker.npy")
  out_inv = outter +1
  out_inv[np.where(out_inv>1)]=0
  inner = np.load("lv_inner_marker.npy")
  out_marker = outter[40:200,100:,150:]
  in_marker = inner[40:200,100:,150:]
#  for i in range (50,60):
#      fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))#figure(figsize=(16,8))
#      axes.imshow(out_marker[i]*2+in_marker[i])
#      #axes[1].imshow(in_marker[i]*2)
#      fig.show() 
  mask = (in_marker + out_marker) +1
  mask[np.where(mask>1)]=0 
  filename = "lv_pot.npy"
  #cc3df.generate_potential(mask, inner[40:200,100:,150:], out_inv[40:200,100:,150:], filename)
  pot = np.load(filename)
  dimb = inner.shape
  final_size = np.zeros(dimb)-1
  final_size[40:200,100:,150:] = pot
  for i in range (6):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    axes.imshow(final_size[40+i*30], cmap = plt.cm.jet)
    #fig.savefig("Potential"+str(40 + i*30)+".png")
    fig.show()
  #np.save("full_lv_potential.npy", final_size)
    
if False :
    from skimage.segmentation import random_walker
    from skimage.morphology import binary_dilation, binary_erosion
    #generate heart potential
    lv = open_mha("C:/Users/cjaquet/Documents/SynthTreeData/3c9e679d-2eab-480c-acaa-31da12301b0a/lv.mha")
    heart = open_mha("C:/Users/cjaquet/Documents/SynthTreeData/3c9e679d-2eab-480c-acaa-31da12301b0a/heart.mha")  
    dim_h = lv.shape

    heart_dil = binary_dilation(heart[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50], np.ones((15,15,15)))    
    heart_mask=heart_dil+1
    heart_mask[np.where(heart_mask>1)]=0     
    lv_mask = binary_erosion(lv[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50], np.ones((5,5,5)))
    mask = (lv_mask + heart_mask) 

    
    inner_marker = lv[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50]
    heart_small_dil =  binary_dilation(heart, np.ones((10,10,10)))
    heart_marker = heart_small_dil +1
    heart_marker[np.where(heart_marker>1)]=0 
    outer_marker = heart_marker[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50]
    filename = "heart_pot.npy"
    #result = cc3df.generate_potential_d(mask, inner_marker, outer_marker, filename)
    markers = outer_marker*2 + inner_marker
    result = random_walker(mask, markers, copy =True, return_full_prob = True, mode= 'cg_mg')
    result[0][np.where(result[0] == 0.)] = -1
    
    final_size = np.zeros(dim_h)-1
    final_size[10:int(dim_h[0]) -2,100:int(dim_h[1])-65, 32:int(dim_h[2])-50] = result[0]
    np.save("heart_potential.npy", final_size)

from skimage.morphology import medial_axis    
if False:
    inner = np.load("lv_inner_marker.npy")
    # Compute the medial axis (skeleton) and the distance transform
    skel, distance = medial_axis(inner, return_distance=True)
    # Distance to the background for pixels of the skeleton
    dist_on_skel = distance * skel
    print dist_on_skel
    
    