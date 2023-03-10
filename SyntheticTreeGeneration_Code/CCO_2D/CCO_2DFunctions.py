# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:18:32 2016

@author: cjaquet
"""

import numpy as np
import sys
import copy

#########################################################

## distance criteria
# n_term is the final number of terminal segment
# k_term is the current number of terminal segment (at the current growth step)

#the distance criterion d_tresh is defined by k_term and r_supp (microbox radius) 
def calculate_d_tresh_2D(r_supp, k_term):
    r_pk = np.sqrt((k_term + 1) * r_supp**2)
    return np.sqrt(np.pi*(r_pk**2)/(k_term)), r_pk
    
## random location generation
def random_location():
    return np.random.rand(2)*200
       
    
# area is defined by two descriptors: [0] = center coord(x,y) and [1] = radius
def belongs_to_area(coords, area):
    vector = coords-area[0]
    if np.sqrt(vector[0]**2 + vector[1]**2) <= area[1] :
        return True
    else:
        return False

# the location of the first segment end is not constrained by the distance criterion, only by the perfusion territory
def first_segmt_end(area, area_descptr):
    inside_area = False
    while inside_area == False :    
        position = random_location()
        if (belongs_to_area(position, area_descptr)):
            return position


#######################################################

def test_connection_list(list_input):
    copy_tree = copy.deepcopy(list_input[0])
    return copy_tree.test_connection(list_input[1], list_input[2])


#######################################################
#radius study
            
def collect_radii_along_bifurcations(tree, filename):
    radii = []
    indexes = [[tree.get_root_index()]]
    a = 0
    all_leaves_reached = False
    while (all_leaves_reached == False):
        children = []
        radius_list = []
        for i in indexes[a]:
            children_array = tree.nodes[i].children_index
            if children_array[0] > 0:
                children.append(children_array[0])
            if children_array[1] > 0:
                children.append(children_array[1])
            radius_list.append(tree.get_radius(i))
        radii.append(radius_list)
        if len(children) > 0:
            indexes.append(children)
            a = a+1
        else:
            all_leaves_reached = True
    np.save(filename,radii)

#######################################################

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
        return np.sqrt(vect_ap[0]**2 + vect_ap[1]**2)
    
    relative_position = (vect_ap[0]*vect_ab[0] + vect_ap[1]*vect_ab[1]) / squared_length_ab
    if (relative_position < 0.): # outside segment passed point a    
        return  np.sqrt(vect_ap[0]**2 + vect_ap[1]**2)
    elif (relative_position > 1.): #outside the segment passed point b     
        return np.sqrt(vect_bp[0]**2 + vect_bp[1]**2)

    projected_point = seg_pt_a + relative_position * vect_ab
    vect_projp = point - projected_point
    return np.sqrt(vect_projp[0]**2 + vect_projp[1]**2)
  
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

#test if segment ab intersect with cd considering their respective width
def no_overlap(point_a, point_b, point_c, point_d, width_ab, width_cd):
    intersect, pos = determinant(point_a, point_b, point_c, point_d)
    if (intersect): 
        #print "intersection between lines"
        if point_is_inside_rectangle(point_a, point_b, width_ab, pos) and point_is_inside_rectangle(point_c, point_d, width_cd, pos):
            #print "intersection is inside segments"
            return False        
        
        dist_a_pos = np.sqrt(np.sum((pos - point_a)**2)) 
        dist_b_pos = np.sqrt(np.sum((pos - point_b)**2)) 
        dist_c_pos = np.sqrt(np.sum((pos - point_c)**2)) 
        dist_d_pos = np.sqrt(np.sum((pos - point_d)**2)) 
        radii_sum = width_ab + width_cd 

        seg_ab_close_to_intersection = False
        seg_cd_close_to_intersection = False
        if dist_a_pos < radii_sum or dist_b_pos < radii_sum:
            seg_ab_close_to_intersection = True
        
        if dist_c_pos < radii_sum or dist_d_pos < radii_sum:
            seg_cd_close_to_intersection = True
            
        if seg_ab_close_to_intersection and seg_cd_close_to_intersection:
            #print "both segment close to intersection"
            return False
            
        # combinaison:
        if point_is_inside_rectangle(point_a, point_b, width_ab, pos) and seg_cd_close_to_intersection:
            #print "intersection inside rectangle and other seg close to intersection"
            return False
        
        if point_is_inside_rectangle(point_c, point_d, width_cd, pos) and seg_ab_close_to_intersection:
            #print "intersection inside rectangle and seg close to intersection"
            return False

    else:
        #print "no intersection between lines"
        #parallel vectors --> no intersection (#vec_ab[0]/vec_cd[0] == vec_ab[1]/vec_cd[1])
        #get distance between segments
        if (segment_distance(point_a, point_b, point_c) < (width_ab + width_cd)):
            return False
        if (segment_distance(point_a, point_b, point_d) < (width_ab + width_cd)):
            return False
    return True
    
def point_is_inside_rectangle(seg_pt_a, seg_pt_b, width, point):   
    #rectangle definition
    vec_ab = seg_pt_b - seg_pt_a 
    
    normal_to_ab = np.array([-vec_ab[1], vec_ab[0]])
    norm_to_ab = normalize(normal_to_ab)
    a = seg_pt_a - norm_to_ab * 0.5 * width
    b = seg_pt_a + norm_to_ab * 0.5 * width
    #c = seg_pt_b + norm_to_ab * 0.5 * width
    d = seg_pt_b - norm_to_ab * 0.5 * width
    
    #mesuring
    m = point
    vec_am = m - a
    vec_ab = b - a
    vec_ad = d - a
    
    ab_am = dot(vec_ab, vec_am)
    ad_am = dot(vec_ad, vec_am)
    ab_ab = dot(vec_ab, vec_ab)
    ad_ad = dot(vec_ad, vec_ad)
    
    if ab_am > 0 and ab_am < ab_ab:
            if ad_am > 0 and ad_am < ad_ad:
                return True
    return False
        
    
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
    


