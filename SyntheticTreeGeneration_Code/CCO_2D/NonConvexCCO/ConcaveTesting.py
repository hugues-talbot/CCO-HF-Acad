# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:46:48 2016

@author: cjaquet
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage.segmentation import random_walker
from scipy.interpolate import RegularGridInterpolator

#### testing random walker to create potential between two surfaces #####

##image definition
if False:
    ext_radius = 200
    int_radius = 50
    a = np.zeros((512, 512)).astype('uint8')
    cx, cy = 256, 256 # The center of circle
    y, x = np.ogrid[-ext_radius: ext_radius, -ext_radius: ext_radius]
    index = x**2 + y**2 <= ext_radius**2
    a[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 1
    
    y_i, x_i = np.ogrid[-int_radius: int_radius, -int_radius: int_radius]
    index_int = x_i**2 + y_i**2 <= int_radius**2
    a[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 0
      
    #markers definition
    markers = np.zeros((512, 512)).astype('uint8')
    index = np.logical_and(x**2 + y**2 >= ext_radius**2 * 0.95, x**2 + y**2 <= ext_radius**2)
    markers[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 2
    index_int = np.logical_and(x_i**2 + y_i**2 >= int_radius**2 * 0.85, x_i**2 + y_i**2 <= int_radius**2)
    markers[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 3

    result = random_walker(a, markers, copy =True, return_full_prob = True)
    
    plt.subplot(1,4, 1)
    plt.imshow(a)
    plt.subplot(1,4,2)
    plt.imshow(markers)
    plt.subplot(1,4,3)
    plt.imshow(result[0])
    plt.subplot(1,4,4)
    plt.imshow(result[1])
    plt.show()
    
def potential_image(center, ext_radius, int_radius):
    cx,cy = center[0], center[1]
    im = np.zeros((cx+ext_radius*2, cy+ext_radius*2)).astype('uint8')
    markers = np.zeros((cx+ext_radius*2, cy+ext_radius*2)).astype('uint8')

    y,x = np.ogrid[-ext_radius: ext_radius, -ext_radius: ext_radius]   
    index = x**2 + y**2 <= ext_radius**2
    print index
    im[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 1
    
    y_i, x_i = np.ogrid[-int_radius: int_radius, -int_radius: int_radius]
    index_int = x_i**2 + y_i**2 <= int_radius**2
    im[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 0
    
    index = np.logical_and(x**2 + y**2 >= ext_radius**2 * 0.94, x**2 + y**2 <= ext_radius**2)
    markers[cy-ext_radius:cy+ext_radius, cx-ext_radius:cx+ext_radius][index] = 2
    index_int = np.logical_and(x_i**2 + y_i**2 >= int_radius**2 * 0.84, x_i**2 + y_i**2 <= int_radius**2)
    markers[cy-int_radius:cy+int_radius, cx-int_radius:cx+int_radius][index_int] = 3
    
    result = random_walker(im, markers, copy =True, return_full_prob = True)
    
    plt.imshow(result[1])
    return result[1]
    

def get_potential(potential, location):
    #find interpolation solution
    x,y = np.arange(0,result.shape[0],1), np.arange(0,result.shape[1],1)
    my_interpolating_function = RegularGridInterpolator((x,y),potential)   
    return my_interpolating_function(location[::-1])
    
def get_interpolating_function(grid):
    x,y = np.arange(0,grid.shape[0],1), np.arange(0,grid.shape[1],1)
    return RegularGridInterpolator((x,y),grid)
    
ext_radius = 50
int_radius = 15
center=np.array([100, 50])
w = potential_image(center, ext_radius, int_radius)
gy,gx = np.gradient(result)

#testing gradient and interpolations
interp_w = get_interpolating_function(w)
interp_gx = get_interpolating_function(gx)
interp_gy = get_interpolating_function(gy)
loc = np.array([100,30])
print "w(loc)", w[loc[1], loc[0]]
print "interpw(loc)", interp_w(loc[::-1])
print "gx(loc)", gx[loc[1], loc[0]]
print "interpgx", interp_gx(loc[::-1])
print "gy", gy[loc[1], loc[0]]
print "interpgy", interp_gy(loc[::-1])

#gdt_vec = np.array([interp_gx(midp[::-1]), interp_gy(midp[::-1])])
def starting_point(interp_w,interp_gx, interp_gy, seg_pt1, seg_pt2, new_location, eps):
    #mid point of seg
    midp = (seg_pt1 + seg_pt2)*0.5
    #find gdt of w at mid point
    gdt_vec = np.array([interp_gx(midp[::-1]), interp_gy(midp[::-1])])
    #find on this line the point p where w(p) = 0.5 * (w(seg_pt1) + w(seg_pt2))
    target_w = 0.5 * (interp_w(seg_pt1[::-1]) + interp_w(seg_pt2[::-1]))
    start_point, w_start = newton_algo(interp_w, interp_gx, interp_gy, new_location, gdt_vec, target_w, eps)
    return start_point

def length(vect):
    return np.sqrt(np.sum(vect**2))

target_line_vect = np.array([0.5,0.5])
target_w = 0.2
point = loc
def newton(w,interp_w, interp_gx, interp_gy, point, target_line_vect, target_w):
    p_gdt = np.array([[interp_gx(point[::-1])], [interp_gy(point[::-1])]])
    #print "p_gdt", p_gdt
    #print p_gdt.shape
    #print "length p_gdt", length(p_gdt)
    #print "interp at point", interp_w(point[::-1])
    norm_line_vect = target_line_vect * 1./ length(target_line_vect)
    
    vect_proj_length = np.sum(p_gdt*target_line_vect) * 1./length(target_line_vect)
    #print "vect_proj length", vect_proj_length
    vect_proj = norm_line_vect * vect_proj_length
    #print "vect_proj", vect_proj
    scal_gap = - interp_w(point[::-1]) + interp_w((point + vect_proj)[::-1])
    #print "scal_gap", scal_gap
    scal_gap_2 = target_w - interp_w(point[::-1])
    #print "scal_gap_2", scal_gap_2
    k = (scal_gap_2 / scal_gap )
    
    scaling = vect_proj_length 
    proj_pt = point + k*scaling * norm_line_vect
    #print "k", k
    #print "norm line vect", norm_line_vect
    print "proj pt",  proj_pt
    #print "w at proj point", interp_w(proj_pt[::-1])
    fac = 10.
    plt.subplot(1,1,1)
    plt.imshow(w)
    plt.plot(point[0],point[1], "o", color = "k", markersize=1)
    plt.plot([point[0], point[0] + p_gdt[0]*fac], [point[1], point[1] + p_gdt[1]*fac], color = "k")   
    plt.plot(proj_pt[0], proj_pt[1], "o", color = "r", markersize = 2)
    plt.plot([point[0], point[0] + target_line_vect[0]*fac], [point[1], point[1] + target_line_vect[1]*fac], color = "b")
    plt.plot([point[0], point[0] + vect_proj[0]*fac], [point[1], point[1] + vect_proj[1]*fac], color = "y")
    plt.show()
    return proj_pt
    
a = newton(w,interp_w, interp_gx, interp_gy, loc, target_line_vect, target_w)

def newton_step(interp_w, interp_gx, interp_gy, point, target_line_vect, target_w):
    p_gdt = np.array([[interp_gx(point[::-1])], [interp_gy(point[::-1])]])
    norm_line_vect = target_line_vect * 1./ length(target_line_vect)  
    vect_proj_length = np.sum(p_gdt*target_line_vect) * 1./length(target_line_vect)   
    vect_proj = norm_line_vect * vect_proj_length
    scal_gap = - interp_w(point[::-1]) + interp_w((point + vect_proj)[::-1])
    scal_gap_2 = target_w - interp_w(point[::-1])
    k = (scal_gap_2 / scal_gap )   
    scaling = vect_proj_length 
    proj_pt = point + k*scaling * norm_line_vect
    return proj_pt
    
def newton_algo(interp_w, interp_gx, interp_gy, point, target_line_vect, target_w, eps):
    proj_pt = newton_step(interp_w, interp_gx, interp_gy, point, target_line_vect, target_w)
    w_proj = interp_w(proj_pt[::-1])
    print "w_proj", w_proj
    if w_proj < target_w + eps and w_proj > target_w - eps:
        return proj_pt, w_proj
    else:
        newton_algo(interp_w, interp_gx, interp_gy, proj_pt, target_line_vect, target_w, eps)

eps = 0.001  
newton_algo(interp_w, interp_gx, interp_gy, loc, target_line_vect, target_w, eps)
#newton(w, interp_w, interp_gx, interp_gy, loc, target_line_vect,target_w)