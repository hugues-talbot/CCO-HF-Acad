# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:18:32 2016

@author: cjaquet
"""

import numpy as np
import sys
import copy
from skimage.segmentation import random_walker
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl
import mpl_toolkits.mplot3d as a3

#from pyamg import *

###########################################Non convex functions####################
## creating a potential grid using random walker
#inner surface and concavity potential value is 1
#outer surface potential value is 0

    
def generate_potential(mask, inner_marker,outer_marker, filename):
    markers = outer_marker*2 + inner_marker
    result = random_walker(mask, markers, copy =True, return_full_prob = True, mode= 'cg_mg')
    result[1][np.where(result[1] == 0.)] = -1
    if (filename != ""):
        np.save(filename,result[1])
    return result[1]
    

def potential_ellipsoid(center, r_int, r_ext, half,cut_top, cutop_val):
    print "starting potential generation"
    #initialization
    rz,ry,rx = int(r_ext[0]), int(r_ext[1]), int(r_ext[2])# r_ext[0] is smallest radius
    riz,riy,rix = int(r_int[0]), int(r_int[1]), int(r_int[2])#r_int[0] is smallest radius
    cx,cy,cz = int(center[0]), int(center[1]),int(center[2])
    im = np.zeros((cx+rx*2, cy+rx*2,cz+rx*2)).astype('uint8')
    markers = np.zeros((cx+rx*2, cy+rx*2, cz+rx*2)).astype('uint8')
    margin = int(np.ceil(rx/10.))
    # grids of index   
    z,y,x = np.ogrid[-rz-margin: rz+margin, -ry-margin: ry+margin, -rx-margin:rx+margin]  
    z_i, y_i, x_i = np.ogrid[-riz: riz, -riy: riy,-rix: rix]
    
    #mask creation
    # a point is inside an ellipsoid if (x/a)2 + (y/b)2 + (z/c)2 <= 1
    index = (x/float(rx))**2 + (y/float(ry))**2 + (z/float(rz))**2 <= 1
    im[cz-rz-margin:cz+rz+margin, cy-ry-margin:cy+ry+margin, cx-rx-margin:cx+rx+margin][index] = 1
    index_int = (x_i/float(rix))**2 + (y_i/float(riy))**2 + (z_i/float(riz))**2 <= 1   
    im[cz-riz:cz+riz, cy-riy:cy+riy, cx-rix:cx+rix][index_int] = 0

    #marker creation
    #external
    index = (x/float(rx))**2 + (y/float(ry))**2 + (z/float(rz))**2 > 1   
    markers[cz-rz-margin:cz+rz+margin, cy-ry-margin:cy+ry+margin, cx-rx-margin:cx+rx+margin][index] = 2
    #internal
    index_int = (x_i/float(rix))**2 + (y_i/float(riy))**2 + (z_i/float(riz))**2 <= 1
    markers[cz-riz:cz+riz, cy-riy:cy+riy, cx-rix:cx+rix][index_int] = 3

    #random_walker    
    print "starting random walker"   
    result = random_walker(im, markers, copy =True, return_full_prob = True, mode= 'cg_mg')   
    print "rdm walker shape",result.shape
    
    #updating final shape: 
    result[1][np.where(result[1] == 0.)] = -1
    if half :
        result[1][0:center[0],:,:] = -1
    if cut_top:
        thresh = center[2]+rx - cutop_val
        result[1][:,:,thresh:result[1].shape[2]] = -1
    
    #printing option
#    slicing = center[2]+5
#    ax = a3.Axes3D(pl.figure(figsize=(8, 8)))
#    xx, yy = pl.ogrid[0:im.shape[1], 0:im.shape[2]]
#    #xx, yy = np.meshgrid(np.linspace(0,1,12), np.linspace(0,1,13))
#    # create vertices for a rotated mesh (3D rotation matrix)
#    X =  xx 
#    Y =  yy
#    Z =  slicing*np.ones(X.shape)
#
#    #cmap for random walker
#    N=50
#    colors = [(1.0,1.0,1.0)]
#    colors.extend(plt.cm.jet(np.linspace(0., 1., N)))
#    colors.extend([(1.0,1.0,1.0)])
#    cmap =mpl.colors.ListedColormap(colors) #plt.cm.jet
#    
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(result[1][slicing,:,:,].transpose()),shade=False)
#    slicing = center[2] +1
#    Z =  slicing*np.ones(X.shape)
#    ax.plot_surface(Z,Y,X, rstride=1, cstride=1, facecolors=cmap(result[1][slicing,:,:,].transpose()),shade=False)
#    ax.set_xlabel('X axis')
#    ax.set_ylabel('Y axis')
#    ax.set_zlabel('Z axis')
#    plt.show()
    
    return result[1]


def length(vect, vox_size):
    if vox_size[0] == 0. and vox_size[1] == 0. and vox_size[2] == 0.:
        return  np.sqrt(np.sum((vect)**2))
    return np.sqrt(np.sum((vect*vox_size)**2))
    
#########################################################

## distance criteria
# n_term is the final number of terminal segment
# k_term is the current number of terminal segment (at the current growth step)

#v_perf is the total perfused volume     
def calculate_d_tresh_3D(r_supp, k_term): 
    r_pk = np.power((k_term + 1)* r_supp**3,3)
    return np.power(np.pi*4.*np.power(r_pk,3)/(3.*k_term),1./3.),r_pk
    
## random location generation
def random_location(im_size):
    position = np.random.rand(3)*im_size
    return position
       

#######################################################

def test_connection_list(list_input):
    #copy_tree = copy.deepcopy(list_input[0])
    return list_input[0].test_connection(list_input[1], list_input[2],list_input[3], list_input[4], list_input[5])


#########################################################

#sibling ratio is r_child_0 / r_child_1
def calculate_betas(sibling_ratio, gamma):
    beta_child_0 = np.power(1 + sibling_ratio**-gamma, -1./gamma)
    beta_child_1 = np.power(1 + sibling_ratio**gamma, -1./gamma)
    #print "beta child_0", beta_child_0, "beta child_1", beta_child_1
    return np.array([beta_child_0, beta_child_1])
    
#return distance between point and segment
def segment_distance(seg_pt_a, seg_pt_b, point):
    vect_ab = seg_pt_b - seg_pt_a
    vect_ap = point - seg_pt_a
    vect_bp = point - seg_pt_b
    squared_length_ab = float(np.sum(vect_ab**2))
    if (squared_length_ab == 0.): 
        return np.sqrt(np.sum(vect_ap**2))
    
    relative_position = (np.sum(vect_ap*vect_ab)) / squared_length_ab
    if (relative_position < 0.): # outside segment passed point a    
        return  np.sqrt(np.sum(vect_ap**2))
    elif (relative_position > 1.): #outside the segment passed point b     
        return np.sqrt(np.sum(vect_bp**2))

    projected_point = seg_pt_a + relative_position * vect_ab
    vect_projp = point - projected_point
    return np.sqrt(np.sum(vect_projp**2))
  
#return distance between line and point  
def line_distance(seg_pt_a, seg_pt_b, point):
    vec_ab = seg_pt_b - seg_pt_a
    squared_length_ab = np.sum(vec_ab**2)
    
    #a and b on top of each other
    if (squared_length_ab < sys.float_info.epsilon):
        return -1.
    
    vec_apoint = point - seg_pt_a
    vec_ab_dot_vec_apoint = np.sum(vec_apoint * vec_ab)
    relative_position = vec_ab_dot_vec_apoint / squared_length_ab
    
    projected_point = seg_pt_a + relative_position * vec_ab
    normal_vec = projected_point - point
    return np.sqrt(np.sum(normal_vec**2))

#test if any of the 3 segments has a length close to 0 (the bifurcation has degenerated into 2 segments)  
def degenerating_test(c0, c1, c2, branching_location, radii, length_factor):
    seg_parent_length = np.sqrt(np.sum((c0 - branching_location)**2)) * length_factor
    seg_old_child_length = np.sqrt(np.sum((c1 - branching_location)**2)) * length_factor
    seg_new_child_length = np.sqrt(np.sum((c2 - branching_location)**2)) * length_factor

    if (seg_parent_length < 2.*radii[0]):
        #print "parent seg degenerate"            
        return False

    if (seg_old_child_length < 2.*radii[1]):
        #print "old child seg degenerate"
        return False
 
    if (seg_new_child_length < 2.*radii[2]):
        #print "new child seg degenerate"
        return False
    
    return True


    
def closestDistanceBetweenLines(a0,a1,b0,b1,clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):
    ''' Given two line segments defined by numpy.array pairs (a0,a1,b0,b1)
        Return the two closest points, and the distance
    '''

    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True

    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0

    _A = A / np.linalg.norm(A)
    _B = B / np.linalg.norm(B)
    cross = np.cross(_A, _B);

    denom = np.linalg.norm(cross)**2


    # If denominator is 0, lines are parallel: Calculate distance with a projection
    # and evaluate clamp edge cases
    if (denom == 0):
        d0 = np.dot(_A,(b0-a0))
        d = np.linalg.norm(((d0*_A)+a0)-b0)

        # If clamping: the only time we'll get closest points will be when lines don't overlap at all.
        # Find if segments overlap using dot products.
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))

            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 == True and clampB1 == True:
                    if np.absolute(d0) < np.absolute(d1):
                        return b0,a0,np.linalg.norm(b0-a0)
                    return b1,a0,np.linalg.norm(b1-a0)

            # Is segment B after A?
            elif d0 >= np.linalg.norm(A) <= d1:
                if clampA1 == True and clampB0 == True:
                    if np.absolute(d0) < np.absolute(d1):
                        return b0,a1,np.linalg.norm(b0-a1)
                    return b1,a1,np.linalg.norm(b1,a1)

        # If clamping is off, or segments overlapped, we have infinite results, just return position.
        return None,None,d


    # Lines criss-cross: Calculate the dereminent and return points
    t = (b0 - a0);
    det0 = np.linalg.det([t, _B, cross])
    det1 = np.linalg.det([t, _A, cross])

    t0 = det0/denom;
    t1 = det1/denom;

    pA = a0 + (_A * t0);
    pB = b0 + (_B * t1);

    # Clamp results to line segments if needed
    if clampA0 or clampA1 or clampB0 or clampB1:

        if t0 < 0 and clampA0:
            pA = a0
        elif t0 > np.linalg.norm(A) and clampA1:
            pA = a1

        if t1 < 0 and clampB0:
            pB = b0
        elif t1 > np.linalg.norm(B) and clampB1:
            pB = b1

    d = np.linalg.norm(pA-pB)

    return pA,pB,d
    
#test if segment ab intersect with cd considering their respective width
def no_overlap(point_a, point_b, point_c, point_d, width_ab, width_cd):
    p1,p0,dist = closestDistanceBetweenLines(point_a, point_b,point_c,point_d, True)
    if dist < (width_ab + width_cd):
        return False
    else :
        return True
        
    
def normalize(vector):
    length = np.sqrt(np.sum(vector**2, axis = 0))
    if length != 0:
        return vector/length
    else:
        return vector
        
def dot(vec_a, vec_b):
    return vec_a[0]*vec_b[0] + vec_a[1]*vec_b[1]    
    
    
def determinant(pt_a, pt_b,pt_c,pt_d):
    denom = (pt_a[0] - pt_b[0])*(pt_c[1] - pt_d[1]) - (pt_a[1] - pt_b[1])*(pt_c[0] - pt_d[0])
    if (denom == 0.):
        return False, np.empty(2)
    x = ( (pt_a[0]*pt_b[1] - pt_a[1]*pt_b[0])*(pt_c[0]-pt_d[0]) - (pt_a[0] - pt_b[0])*(pt_c[0]*pt_d[1] - pt_c[1]*pt_d[0]) ) / denom 
    y = ( (pt_a[0]*pt_b[1] - pt_a[1]*pt_b[0])*(pt_c[1]-pt_d[1]) - (pt_a[1] - pt_b[1])*(pt_c[0]*pt_d[1] - pt_c[1]*pt_d[0]) ) / denom
    return True, np.array([x,y])
    


