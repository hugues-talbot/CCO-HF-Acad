# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:18:32 2016

@author: cjaquet
"""

import numpy as np
import sys
import copy
from skimage.segmentation import random_walker



#from pyamg import *

###########################################Non convex functions####################
## creating a potential grid using random walker
#inner surface and concavity potential value is 1
#outer surface potential value is 0

    
def generate_potential(mask, inner_marker,outer_marker):
    markers = outer_marker*2 + inner_marker
    result = random_walker(mask, markers, copy =True, return_full_prob = True, mode= 'cg_mg')
    result[1][np.where(result[1] == 0.)] = -1
#    if (filename != ""):
#        np.save(filename,result[1])
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
def segment_distance(seg_pt_a, seg_pt_b, point, vox_size):
    vect_ab = seg_pt_b - seg_pt_a
    vect_ap = point - seg_pt_a
    vect_bp = point - seg_pt_b
    squared_length_ab = length(vect_ab, vox_size)
    if (squared_length_ab == 0.): 
        return length(vect_ap, vox_size)
    
    relative_position = (np.sum(vect_ap*vect_ab)) / squared_length_ab
    if (relative_position < 0.): # outside segment passed point a    
        return  length(vect_ap, vox_size)
    elif (relative_position > 1.): #outside the segment passed point b     
        return length(vect_bp, vox_size)

    projected_point = seg_pt_a + relative_position * vect_ab
    vect_projp = point - projected_point
    return length(vect_projp, vox_size)
  
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
def degenerating_test(c0, c1, c2, branching_location, radii, length_factor, vox_size):
    seg_parent_length = length(c0 - branching_location, vox_size) * length_factor
    seg_old_child_length = length(c1 - branching_location, vox_size) * length_factor
    seg_new_child_length = length(c2 - branching_location, vox_size) * length_factor
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
            #print "case 1"
            pA = a0
        elif t0 > np.linalg.norm(A) and clampA1:
            #print "case 2"
            pA = a1

        if t1 < 0 and clampB0:
            #print "case 3"
            pB = b0
        elif t1 > np.linalg.norm(B) and clampB1:
            #print "case 4"
            pB = b1

    d = np.linalg.norm(pA-pB)

    return pA,pB,d
    
#test if segment ab intersect with cd considering their respective width
def no_overlap(point_a, point_b, point_c, point_d, width_ab, width_cd, vox_size):
    p1,p0,dist = closestDistanceBetweenLines(point_a, point_b,point_c,point_d, True)
        
    #if np.all(dist / vox_size > (width_ab + width_cd)/vox_size):
    if dist > (width_ab + width_cd):
        #print "dist / vox size", dist / vox_size, "width (width_ab + width_cd)/vox_size", (width_ab + width_cd)/vox_size
        print "p1", p1 , "p0", p0 , "dist ", dist , "(width_ab + width_cd)", (width_ab + width_cd)
        return True
    else:
        return False
#    if dist < (width_ab + width_cd):
#        return False
#    else :
#        return True
        
    
def normalize(vector):
    length = np.sqrt(np.sum(vector**2, axis = 0))
    if length != 0:
        return vector/length
    else:
        return vector
        
def dot(vec_a, vec_b):
    return np.sum(vec_a*vec_b)    
    
    
def determinant(pt_a, pt_b,pt_c,pt_d):
    denom = (pt_a[0] - pt_b[0])*(pt_c[1] - pt_d[1]) - (pt_a[1] - pt_b[1])*(pt_c[0] - pt_d[0])
    if (denom == 0.):
        return False, np.empty(2)
    x = ( (pt_a[0]*pt_b[1] - pt_a[1]*pt_b[0])*(pt_c[0]-pt_d[0]) - (pt_a[0] - pt_b[0])*(pt_c[0]*pt_d[1] - pt_c[1]*pt_d[0]) ) / denom 
    y = ( (pt_a[0]*pt_b[1] - pt_a[1]*pt_b[0])*(pt_c[1]-pt_d[1]) - (pt_a[1] - pt_b[1])*(pt_c[0]*pt_d[1] - pt_c[1]*pt_d[0]) ) / denom
    return True, np.array([x,y])
    


