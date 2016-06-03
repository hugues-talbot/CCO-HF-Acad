# -*- coding: utf-8 -*-
"""
Created on Tue May 24 11:26:13 2016

@author: cjaquet
"""
import numpy as np
import sys


#########################################################

## distance criteria
# n_term is the final number of terminal segment
# k_term is the current number of terminal segment (atg the current growth step)

#v_perf is the total perfused volume    
def calculate_r_supp_3D(v_perf, n_term):
    return np.power(v_perf*3./(n_term*4.*np.pi), 1./3.)
     
def calculate_d_tresh_3D(r_supp, k_term):    
    return np.power(np.pi*4.*np.power(r_supp,3)/(3.*k_term),1./3.)
    
def get_d_tresh(v_perf, n_term, k_term):
    return calculate_d_tresh_3D(calculate_r_supp_3D(v_perf, n_term), k_term)
    
## random location generation
def random_location():
    return np.random.rand(3)*200
    
## test if new location is over d_tresh distance from existing segments
def test_dist_criteria(tree, location, d_tresh, area):
    if (belongs_to_volume(location, area)):    
        for sgmt in tree.nodes:
            if (sgmt.parent() >= 0):
                dist = segment_distance(sgmt.coord, (tree.get_node(sgmt.parent())).coord, location)
                if (dist < d_tresh):
                    #print "point too close of existing segments"
                    #print "dist", dist, "d_tresh", d_tresh
                    #print "location", location
                    return False
    else:
        #print "point out of area"
        #print "location", location
        return False
    print "point inside area and over distance threshold"
    return True  
    
    
# area is defined by two descriptors: [0] = center coord(x,y) and [1] = radius
def belongs_to_area(coords, area):
    vector = coords-area[0]
    if np.sqrt(vector[0]**2 + vector[1]**2) <= area[1] :
        return True
    else:
        return False

#volume is defined by two descriptors: [0] = center of sphere coord (x,y,z) and [1] = radius
def belongs_to_volume(coords, volume):
    vector = coords-volume[0]
    if np.sqrt(np.sum(vector**2)) <= volume[1] :
        return True
    else:
        return False


def first_segmt_end(volume):
    inside_area = False
    while inside_area == False :    
        position = random_location()
        #position[2]=50
        if (belongs_to_volume(position, volume)):
            return position
            
def get_new_location(tree, volume_descpt, n_term):   
    v_perf = (volume_descpt[1]**3) * np.pi*4./3
    k_term = tree.get_k_term()
    print "k_term", k_term
    d_tresh = get_d_tresh(v_perf, n_term, k_term)
    print "d_tresh", d_tresh
    meet_criteria = False
    ind = 0
    while (meet_criteria == False and ind < 1000):
        point = random_location()
        #point[2]=50
        if (test_dist_criteria(tree, point, d_tresh, volume_descpt)):
            print "location found"
            return True, point, d_tresh
        ind = ind + 1
        if ind == 1000:
            print "impossible to find new location with current d_tresh"
            d_tresh = 0.9*d_tresh
            print "using new value: ", d_tresh
            ind = 0
    return False, np.array([0.,0.]), d_tresh
    
    
#########################################################

#sibling ratio is r_child_0 / r_child_1
def calculate_betas(sibling_ratio, gamma):
    beta_child_0 = np.power(1 + sibling_ratio**-gamma, -1./gamma)
    beta_child_1 = np.power(1 + sibling_ratio**gamma, -1./gamma)
    print "beta child_0", beta_child_0, "beta child_1", beta_child_1
    return np.array([beta_child_0, beta_child_1])
    
#return distance between point and segment
def segment_distance(seg_pt_a, seg_pt_b, point):
    vect_ab = seg_pt_b - seg_pt_a
    vect_ap = point - seg_pt_a
    vect_bp = point - seg_pt_b
    squared_length_ab = np.sum(vect_ab**2)
    if (squared_length_ab == 0.): 
        print 1
        return np.sqrt(np.sum(vect_ap**2))
    
    relative_position = (np.sum(vect_ap*vect_ab)) / squared_length_ab
    if (relative_position < 0.): # outside segment passed point a    
        return  np.sqrt(np.sum(vect_ap**2))
    elif (relative_position > 1.): #outside the segment passed point b     
        return np.sqrt(np.sum(vect_bp**2))

    projected_point = seg_pt_a + relative_position * vect_ab
    vect_projp = point - projected_point
    return np.sqrt(np.sum(vect_projp**2))
    
def cross(vec_a, vec_b):
    res = np.zeros(3)
    res[0] = (vec_a[1] * vec_b[2]) - (vec_a[2] * vec_b[1])
    res[1] = (vec_a[2] * vec_b[0]) - (vec_a[0] * vec_b[2])
    res[2] = (vec_a[0] * vec_b[1]) - (vec_a[1] * vec_b[0])
    return res

def float_comparison(vec_a, vec_b):
    diff = np.abs(vec_a - vec_b)
    if (diff[0] < sys.float_info.epsilon) and  (diff[1] < sys.float_info.epsilon) and (diff[2] < sys.float_info.epsilon):
        return True
    return False



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
def degenerating_test(c0, c1, c2, branching_location, radii):
    seg_parent_length = np.sum((c0 - branching_location)**2)
    seg_old_child_length = np.sum((c1 - branching_location)**2)
    seg_new_child_length = np.sum((c2 - branching_location)**2)
    #print "branching location", branching_location
    #print "parent length", seg_parent_length, "2*radius", 2.*radii[0]
    if (seg_parent_length < 2.*radii[0]):
        print "parent seg degenerate"            
        return False

    #print "old child length", seg_old_child_length, "2*radius", 2.*radii[1]    
    if (seg_old_child_length < 2.*radii[1]):
        print "old child seg degenerate"
        return False
    
    #print "new child length", seg_new_child_length, "2*radius", 2.*radii[2]   
    if (seg_new_child_length < 2.*radii[2]):
        print "new child seg degenerate"
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
        

        