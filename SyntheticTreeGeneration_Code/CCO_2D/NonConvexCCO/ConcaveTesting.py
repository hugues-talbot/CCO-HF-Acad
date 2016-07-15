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
    


#######testing gradient and interpolations

ext_radius = 50
int_radius = 15
center=np.array([100, 50])
w = potential_image(center, ext_radius, int_radius)
gy,gx = np.gradient(w)

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


def length(vect):
    return np.sqrt(np.sum(vect**2))

def newton(w,interp_w, interp_gx, interp_gy, point, target_line_vect, target_w):
    p_gdt = np.array([interpol(interp_gx,point), interpol(interp_gy,point)])
    print "p_gdt", p_gdt
    print p_gdt.shape
    #print "length p_gdt", length(p_gdt)
    #print "interp at point", interp_w(point[::-1])
    norm_line_vect = target_line_vect * 1./ length(target_line_vect)
    print "norm_line_vect", norm_line_vect
    vect_proj_length = np.sum(p_gdt*target_line_vect) * 1./length(target_line_vect)
    print "vect_proj length", vect_proj_length
    vect_proj = norm_line_vect * vect_proj_length
    print "vect_proj", vect_proj
    scal_gap = - interpol(interp_w,point) + interpol(interp_w,(point + vect_proj))
    print "point + vect proj", point + vect_proj
    print "interp_w((point + vect_proj)[::-1])",interp_w((point + vect_proj)[::-1])
    print "scal_gap", scal_gap
    scal_gap_2 = target_w - interpol(interp_w,point)
    print "scal_gap_2", scal_gap_2
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
    
def newton_algo_visu(w,interp_w, interp_gx, interp_gy, point, target_line_vect, target_w, eps):
    proj_pt = newton(w,interp_w, interp_gx, interp_gy, point, target_line_vect, target_w)
    w_proj = interpol(interp_w,proj_pt)
    print "w_proj", w_proj
    print "up boundary", target_w + eps
    print "down boundary", target_w - eps
    if (w_proj < target_w + eps) and( w_proj > target_w - eps):
        return proj_pt #, w_proj     
    else:
        return newton_algo_visu(w,interp_w, interp_gx, interp_gy, proj_pt, target_line_vect, target_w,eps)

def interpol(interp_func, location):
    return float(interp_func(location[::-1]))
  
def starting_point(w,interp_w,interp_gx, interp_gy, seg_pt1, seg_pt2, new_location, eps):
    #mid point of seg
    midp = (seg_pt1 + seg_pt2)*0.5
    #find gdt of w at mid point
    gdt_vec = np.array([interpol(interp_gx, midp), interpol(interp_gy,midp)])
    #find on this line the point p where w(p) = 0.5 * (w(seg_pt1) + w(seg_pt2))
    target_w = 0.5 * (interp_w(seg_pt1[::-1]) + interp_w(seg_pt2[::-1]))
    print "target w", target_w
    start_point = newton_algo_visu(w,interp_w, interp_gx, interp_gy, new_location, gdt_vec, target_w, eps)
    print "interp starting point", interp_w(start_point[::-1])
    print "start_point", start_point
    #print "w_start", w_start
    plt.imshow(w)
    fac = 10.
    plt.plot([seg_pt1[0], seg_pt2[0]], [seg_pt1[1], seg_pt2[1]], color = "g")
    plt.plot([midp[0], midp[0]+gdt_vec[0]*fac], [midp[1], midp[1]+gdt_vec[1]*fac], color = "b")
    plt.plot(start_point[0], start_point[1], color = "r", markersize = 2)
    plt.plot(new_location[0], new_location[1], color = "k", markersize = 1)
    plt.show()
    return start_point
        
####### testing newton-raphson implementation to find bifurcation point to initialize Kamiya algo with
if False: 
    #target_line_vect = np.array([0.5,0.5])
    target_w = 0.2
    point = loc    
    eps = 0.001  
    #newton_algo_visu(w,interp_w, interp_gx, interp_gy, loc, target_line_vect, target_w, eps)
    seg_1=np.array([80,20])
    seg_2 = np.array([95,20])
    starting_point(w, interp_w, interp_gx, interp_gy, seg_1, seg_2, loc, 0.001)
    
########testing research of sampling n
def local_n(w, interp_w, interp_gx, interp_gy, seg_1, seg_2):
    wp1 = interpol(interp_w,seg_1)
    wp2 = interpol(interp_w,seg_2)
    gdt_p1 = np.array([interpol(interp_gx, seg_1), interpol(interp_gy, seg_1)])
    gdt_p2 = np.array([interpol(interp_gx, seg_1), interpol(interp_gy, seg_2)])
    p1p2_vec = seg_2 - seg_1
    l1 = - wp1 * (1. / np.sum(gdt_p1* p1p2_vec))
    l2 = - wp2 * (1. / np.sum(gdt_p2* -p1p2_vec))
    print "l1", l1, "l2", l2
    lbda = min(l1,l2)
    print "lbda", lbda    
    splg_n = np.ceil(np.abs(1./lbda))
    print "splgn", splg_n
    
    #fac = 10.
    plt.imshow(w)
    norm_vec= p1p2_vec / length(p1p2_vec)
    sample = length(p1p2_vec) / splg_n
    print "sample", sample
    plt.plot([seg_1[0], seg_2[0]], [seg_1[1], seg_2[1]], color = "y")
    locs = []
    for i in range (int(splg_n+1)):
        loc = seg_1 + i * sample * norm_vec
        
        print "loc", loc
        locs.append(loc)
    print locs
    locsa = np.array(locs)
    plt.scatter(locsa[:,0], locsa[:,1], color = "k")
    plt.show()  
    return splg_n
    
def calculate_sampling(tolerance, max_curv_radius, w, interp_w, interp_gx, interp_gy, seg_1, seg_2):
    loc_n = local_n(w, interp_w, interp_gx, interp_gy, seg_1, seg_2)
    r_star = max_curv_radius - tolerance
    print "r_star",r_star
    c= np.sqrt(max_curv_radius**2 - (max_curv_radius-tolerance)**2)
    print c
    global_n = np.ceil(length(seg_2-seg_1) / c)
    print "global n", global_n
    if (loc_n >= global_n):
        print "n is local one", loc_n
        return loc_n
    else:
        locs =[]
        sample = length(seg_1-seg_2) / global_n
        for i in range (int(global_n+1)):
            loc = seg_1 + i * sample * ((seg_2-seg_1)/length(seg_2-seg_1))
        
            print "loc", loc
            locs.append(loc)
   
        locsa = np.array(locs)
        plt.scatter(locsa[:,0], locsa[:,1], color = "r")
        print "final n", global_n
        return global_n

seg_1 = np.array([100,15])
seg_2 = np.array([105,20])    
#find_sampling(w, interp_w, interp_gx, interp_gy, seg_1, seg_2)
calculate_sampling(1,int_radius, w, interp_w, interp_gx, interp_gy, seg_1, seg_2)    