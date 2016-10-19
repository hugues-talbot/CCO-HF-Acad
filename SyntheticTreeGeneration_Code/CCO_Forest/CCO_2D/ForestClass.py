# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:48:47 2016

@author: jaquetc
"""

import numpy as np
import CCO_2DFunctions as cco_2df

import NodeClass as nclass


class Forest:
    def __init__(self, trees, n_term, q_perf, p_drop, visc, a_perf, r_f, w, max_curv_radius,center,real_radius, gamma):
        self.trees = trees
        self.n_term = n_term
        self.final_q_perf = q_perf
        self.q_term = q_perf / float(n_term)
        self.p_drop = p_drop
        self.node_nber = 0 #total of node in forest
        self.nu = visc
        self.a_perf = a_perf
        self.r_supp = np.sqrt(a_perf / (n_term * np.pi)) #microbox radius: average radius of terminal segment perfusion territory when reaching final tree growth
        self.final_perf_radius = r_f #real size of the estimated plain perfusion territory radius (excluding concavity surface)
        self.length_factor = 1 # not a const, is updated during tree growth after each added bifurcation       
        self.max_curv_rad = max_curv_radius
        self.center = center
        self.real_final_radius = real_radius #real size of the final cercle (includes the concavity)
        self.gamma = gamma
        self.w_pot = w

        
    def create_tree(self, nterm, source_location, q_perf):
        tree_index = len(self.trees)
        tree = nclass.Tree(tree_index, [], nterm, q_perf, self.p_drop, self.nu, self.a_perf, self.final_perf_radius, self.w_pot, self.max_curv_rad, self.center,self.real_final_radius, self.gamma)        
        root_node = nclass.Node(0,source_location, q_perf, -1)
        root_node.set_child_0_index(1)
        tree.add_node(root_node)   
        return tree
        
    def add_tree(self, tree):
        self.trees.append(tree)
        self.node_nber = self.node_nber + len(tree.nodes)
        
    def first_segment_end(self):
        first_node_position = (self.trees[0]).first_segmt_end()
        first_node = nclass.Node(1, first_node_position, (self.trees[0]).final_q_perf,0)
        (self.trees[0]).add_node(first_node)
        (self.trees[0]).update_flow()
        self.update_forest_length_factor()
        (self.trees[0]).depthfirst_resistances(0)
        (self.trees[0]).update_volume()
        self.node_nber = self.node_nber +1        
        
        location_found = False
        while location_found == False:          
            result = self.get_new_location(1.)
            location_found = result[0]
        secund_node = nclass.Node(1, result[1], (self.trees[1]).final_q_perf, 0)
        (self.trees[1]).add_node(secund_node)
        (self.trees[1]).update_flow()
        self.update_forest_length_factor()
        self.update_forest_tree_resistances()
        self.update_forest_tree_volumes()
        self.node_nber = self.node_nber +1
        
#    tree.depthfirst_resistances(0)        
#        for tree in self.trees:
#             end = tree.first_segmt_end
#             for other_tree in self.trees:
#                 if other_tree.tree_index != tree.tree_index:
#                     if other_tree
#                     if (cco_2df.no_overlap(i.coord, parent_i.coord, old_child_location, branching_location, radius_i_rescaled, new_branches_radii_rescaled[1]) ==  False):
        
    def get_fk_term(self):
        fk_term = 0
        for tree in self.trees:
            fk_term = fk_term + tree.get_k_term()
        return fk_term
       
    def get_new_location(self, d_tresh_factor):   
        k_term = self.get_fk_term()
        d_tresh, r_pk = cco_2df.calculate_d_tresh_2D(self.r_supp, k_term)
        length_factor = r_pk /  self.final_perf_radius
        d_tresh_factorized = d_tresh / length_factor * d_tresh_factor#d_tresh is calculated in the current k_world dimensions
        meet_criteria = False
        ind = 0
        max_it = 100
        while (meet_criteria == False and ind < max_it):
            point = cco_2df.random_location(self.center, self.real_final_radius)
            respect_criteria = True
            for tree in self.trees:
                if tree.test_dist_criteria(point, d_tresh_factorized) == False:
                    respect_criteria = False
                    break
            if (respect_criteria == True):
                print "location found"
                print point
                return True, point, d_tresh_factorized
            ind = ind + 1
            if ind == max_it:
                d_tresh_factorized = 0.9*d_tresh_factorized
                print "using new value: ", d_tresh_factorized
                ind = 0
        return False, np.array([0.,0.]), d_tresh_factorized
        
    def find_forest_neighbors(self, location,n):
        distances=[]
        dist_type = [("distance", float), ("tree_index", int), ("node_index", int)]
        for j in range(len(self.trees)):
            distances = distances + (self.trees[j]).find_neighbors(location, n)          
        threshold = n if (self.node_nber - len(self.trees)> n) else len(distances)
        print "node_nber", self.node_nber, "Ncon ", n, "len(distances) ", len(distances)
        print "threshold", threshold
        d_array = np.array(distances, dtype = dist_type)
        d_sorted = np.sort(d_array, order = "distance")
        #print "d_sorted",d_sorted
        #print "threshold",threshold
        return [[i[1],i[2]] for i in d_sorted[0 : threshold]] 
        
    #the length factor is a scaling factor for the segment length
    #from the final coordinates world (dimensions of the final perfusion territory)
    #to the current k world (whose perfusion territory size is increased before each added bifurcation)
    def update_forest_length_factor(self):       
        r_pk = np.sqrt((self.get_fk_term() + 1)* self.r_supp**2) #correspond to the plain perfusion territory (excluding concavity surface)
        self.length_factor =  r_pk / self.final_perf_radius
        for tree in self.trees:
            tree.length_factor = self.length_factor 
        print "updating length_factor", self.length_factor   
        
    def update_forest_tree_resistances(self):
        for tree in self.trees:
            tree.depthfirst_resistances(0)
            
    def update_forest_tree_volumes(self):
        for tree in self.trees:
            tree.update_volume()
        
    def test_forest_connection(self, tree_ind, node_index, new_location, debug):       
        current_tree = self.trees[tree_ind]
        print "testing neighbor of tree index", tree_ind
        reslt = current_tree.test_connection(node_index, new_location)
        #1, True, tree_vol, result, old_child_index, new_radii
        current_tree_volume =reslt[2]
        beta_and_bif_location = reslt[3]
        old_child_index=reslt[4]
        new_radii_rescaled = reslt[5]
        old_child_node = current_tree.nodes[old_child_index]
        old_child_coord = old_child_node.coord
        parent_coord = (current_tree.nodes[old_child_node.parent()]).coord
        #need : volume and inputs for intersection test
        if debug:
            print "reslt[0]",reslt[0]
            print "reslt[1]",reslt[1]
            print "parent_coord", parent_coord
            print "current_tree_vbolume", current_tree_volume
            print "new_radii_rescaled",new_radii_rescaled
            print "beta_and_bif_location", beta_and_bif_location
            #return False, 0.,[np.zeros(2), np.zeros(2)], tree_ind, node_index, current_tree_volume
        if reslt[1]== True:
            no_intersection = True
            for tree in self.trees:
                if tree.tree_index != tree_ind:
                    intersect = tree.check_intersection_external_seg(old_child_coord, parent_coord, new_location, beta_and_bif_location[1], new_radii_rescaled)
                    if (intersect == False):
                        no_intersection = False
                        break
            if (no_intersection):
                total_volume = current_tree_volume
                print "current tree volume"
                for tree in self.trees:
                    if tree.tree_index != tree_ind:
                        total_volume = total_volume + tree.tree_volume 
                        print "adding volume of tree", tree.tree_index, " ", tree.tree_volume   
                print "test forest connection: positiv result"
                return True, total_volume, beta_and_bif_location, tree_ind, node_index, current_tree_volume 
            else:
                return False, 0.,[np.zeros(2), np.zeros(2)], tree_ind, node_index, current_tree_volume
        else:
            return False, 0.,[np.zeros(2), np.zeros(2)], tree_ind, node_index, current_tree_volume
    
    def add_connection(self, beta_and_bif_location, tree_index, node_index, new_location,tree_volume):
        if (self.trees[tree_index]).add_connection(node_index, new_location, beta_and_bif_location[1], beta_and_bif_location[0], tree_volume):
            ##update perfusion territory: add a microcirculatory black box
            # by updating the length factor 
            self.update_forest_length_factor()   
            # and updating the resistance on whole tree so all radii are rescaled
            self.update_forest_tree_resistances() 
            self.update_forest_tree_volumes()
            self.node_nber = self.node_nber +2
            return True
        else:
            return False
            
    def printing(self):
        for tree in self.trees:
            print "tree index", tree.tree_index
            tree.printing_full()