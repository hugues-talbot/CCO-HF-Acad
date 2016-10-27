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

        
    def create_tree(self, source_location, q_perf):
        tree_index = len(self.trees)
        nterm = np.ceil(q_perf / self.q_term) #tree_nterm#
        a = (9.6e3-7.98e3)/(250-4000) #0.5* #(pterm1-pterm2)/(nterm1-nterm2) with Karch values: pterm =9.60e3 for nterm = 250  (if working in mmhg then a =-0.0032)
        b = 9.6e3 - a * 250 #72.8
        p_term = a * nterm + b
        print "pterm of tree", p_term
        delta_p = 1.33e4 - p_term
        print "delta p of tree", delta_p
        print "tree of index", tree_index, "nterm",nterm
        print "flow ratio", q_perf/self.final_q_perf
        tree = nclass.Tree(tree_index, [], self.q_term, q_perf, delta_p, self.nu, self.a_perf, self.final_perf_radius, self.w_pot, self.max_curv_rad, self.center,self.real_final_radius, self.gamma)        
        root_node = nclass.Node(0,source_location, q_perf, -1)
        root_node.set_child_0_index(1)
        tree.add_node(root_node)   
        return tree
        
    def add_tree(self, tree):
        self.trees.append(tree)
        self.node_nber = self.node_nber + len(tree.nodes)
        
    def first_segment_end(self):
        dist_max = self.find_sources_dist_max()
        print "dist_max", dist_max
        
        #if different flows need to add first segment on smallest flow tree (first_tree_index need to be found)
        #ifrst forest segment should not be fully random --> need to set maximum distance induced by other existing sources location and flow relationships
        indexes = []
        first_added = 0   
        
        indexes.append(first_added)
        self.trees[first_added].set_activity(True)
        first_node_position = self.first_source_end(self.trees[first_added], dist_max[first_added][0]*dist_max[first_added][1])#(self.trees[first_added]).first_segmt_end()
        #found, location, d_tresh = self.get_new_loc_for_sources(1.,)

        #first_node_position = location  
        first_node = nclass.Node(1, first_node_position, (self.trees[0]).final_q_perf,0)
        (self.trees[first_added]).add_node(first_node)
        (self.trees[first_added]).update_flow()
        self.update_forest_length_factor()
        (self.trees[first_added]).depthfirst_resistances(0)
        (self.trees[first_added]).update_volume()
        self.node_nber = self.node_nber +1        
        
#        for tree in self.trees: #actually should better iterate along decreasing tree flow (because d_tresh will decrease with each added segment)
#            current_index = tree.tree_index            
#            if current_index != first_added:            
#                location_found = False
#                while location_found == False:          
#                    result = self.get_new_location(1.,True) #1 is dtreshfactor
#                    location_found = result[0]
#                secund_node = nclass.Node(1, result[1], (self.trees[current_index]).final_q_perf, 0)
#                (self.trees[current_index]).add_node(secund_node)
#                (self.trees[current_index]).update_flow()
#                self.update_forest_length_factor()
#                self.update_forest_tree_resistances(indexes, True)
#                self.update_forest_tree_volumes(indexes, True)
#                self.node_nber = self.node_nber +1
#                indexes.append(current_index)               
                
        for tree in self.trees:
            current_index = tree.tree_index            
            if current_index != first_added:  
                location_found = False
                while location_found == False:          
                    result = self.get_new_loc_for_sources( 1., tree, dist_max[current_index][0]*dist_max[current_index][1]) #1 is dtreshfactor
                    location_found = result[0]
                secund_node = nclass.Node(1, result[1], (self.trees[current_index]).final_q_perf, 0)               
                (self.trees[current_index]).add_node(secund_node)
                (self.trees[current_index]).update_flow()
                self.update_forest_length_factor()
                self.update_forest_tree_resistances(indexes, True)
                self.update_forest_tree_volumes(indexes, True)
                self.node_nber = self.node_nber +1
                indexes.append(current_index)               
                
        
#    tree.depthfirst_resistances(0)        
#        for tree in self.trees:
#             end = tree.first_segmt_end
#             for other_tree in self.trees:
#                 if other_tree.tree_index != tree.tree_index:
#                     if other_tree
#                     if (cco_2df.no_overlap(i.coord, parent_i.coord, old_child_location, branching_location, radius_i_rescaled, new_branches_radii_rescaled[1]) ==  False):
     
    def new_pos_inside_sqr(self,source_loc, dist_max):
            inf_dist_max = False
            new_pos =np.zeros(2)
            while (inf_dist_max == False):
                rdm_in_sqre = np.random.rand(2)*2.*dist_max
                new_pos = source_loc - (np.ones(2)*dist_max) + rdm_in_sqre
                if cco_2df.length(new_pos-source_loc) < dist_max:  
                    inf_dist_max = True
            return new_pos
        
    def first_source_end(self,tree, max_dist):
        # the location of the first segment end is constrained by perfusion territory and other sources locations
        source_loc = tree.nodes[0].coord        
        inside_area = False
        while inside_area == False :    
            position = self.new_pos_inside_sqr(source_loc, max_dist)
            if (tree.inside_perf_territory(position)):
                n = tree.calculate_sampling(self.max_curv_rad, position, source_loc)
                if tree.sample_and_test(position, source_loc, n) == True:
                    print "first source end found"
                    return position
#                else:
#                    print "first source, location crosses concavity"
#            else:
#                print "first source", "out terrritory"
     
    def get_fk_term(self):
        fk_term = 0
        for tree in self.trees:
            fk_term = fk_term + tree.get_k_term()
        return fk_term
        
    def find_sources_dist_max(self):
        tree_nb = len(self.trees)
        dist_max = np.ones((tree_nb,2))*self.real_final_radius*10.
        #print "dist_max_initial", dist_max
        for tree in self.trees:
            for i in range (tree.tree_index+1, tree_nb):
                dist = cco_2df.length(tree.nodes[0].coord - self.trees[i].nodes[0].coord)
                if dist_max[tree.tree_index][0] > dist:
                    dist_max[tree.tree_index][0] = dist
                    tree_q_perf = tree.final_q_perf
                    i_q_perf = self.trees[i].final_q_perf
                    dist_max[tree.tree_index][1] = tree_q_perf / (i_q_perf + tree_q_perf) 
                if dist_max[i][0] > dist:
                    dist_max[i][0] = dist
                    tree_q_perf = tree.final_q_perf
                    i_q_perf = self.trees[i].final_q_perf
                    dist_max[i][1] = tree_q_perf / (i_q_perf + tree_q_perf)  
        return dist_max
    
    def get_new_loc_for_sources(self, d_tresh_factor, current_tree, dist_max): 
         #d_tresh is used to constrain distance to existing segments while dist_max is used for distance to other sources
        source_coord = current_tree.nodes[0].coord
        k_term = self.get_fk_term()
        d_tresh, r_pk = cco_2df.calculate_d_tresh_2D(self.r_supp, k_term)
        length_factor = r_pk / self.final_perf_radius
        d_tresh_factorized = d_tresh / length_factor * d_tresh_factor#d_tresh is calculated in the current k_world dimensions
        meet_criteria = False
        ind = 0
        max_it = 100
        print "looking for a position respecting dist_max", dist_max, "and d_tresh", d_tresh_factorized
        while (meet_criteria == False and ind < max_it):        
            new_pos = self.new_pos_inside_sqr(source_coord, dist_max)
            if current_tree.inside_perf_territory(new_pos): 
               #print "pos in dist_max criteria", new_pos, "dist_max", dist_max
                respect_criteria = True
                for tree in self.trees:
                    if tree.tree_index != current_tree.tree_index:
                        if tree.test_dist_criteria(new_pos, d_tresh_factorized,True) == False:
                            respect_criteria = False
                            break
                if (respect_criteria == True):
                    #print "new_pos found"
                    #print new_pos
                    n = tree.calculate_sampling(self.max_curv_rad, new_pos, source_coord)
                    if tree.sample_and_test(new_pos, source_coord, n) == True:
                        return True, new_pos, d_tresh_factorized
#                    else:
#                        print "crossing concavity"
            ind = ind + 1
            if ind == max_it:
                d_tresh_factorized = 0.9*d_tresh_factorized
                print "using new value for d_tresh: ", d_tresh_factorized
                ind = 0
           
        return False, np.array([0.,0.]), d_tresh_factorized
    
    def inside_perf_territory(self, location):
        if location[1] < 0. or location[0] < 0.:
            return False
        if location[1] < (self.w_pot.shape[0]-1) and location[0] < (self.w_pot.shape[1]-1): 
            pot_val = round(self.trees[0].get_nearest_w(location),3)
            if (pot_val> 0.) and (pot_val < 1.0): 
                return True
        return False    
       
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
            if (self.inside_perf_territory(point) == True):
                respect_criteria = True
                for tree in self.trees:
                    if tree.test_dist_criteria(point, d_tresh_factorized,False) == False:
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
        print "traditional r_pk",r_pk
        sum_microboxes = self.r_supp**2 #use forest r_supp: average value to get homogeneous distribution        
#        for tree in self.trees:
#            sum_microboxes = sum_microboxes + tree.get_k_term() * tree.r_supp            
#        new_r_pk = np.sqrt(sum_microboxes)
#        print "new r_pk",new_r_pk
        self.length_factor =  r_pk / self.final_perf_radius
        for tree in self.trees:
            tree.length_factor = self.length_factor 
        print "updating length_factor", self.length_factor   
        
    def update_forest_tree_resistances(self,indexes,subset):
        if subset:
            for i in indexes:      
                (self.trees[i]).depthfirst_resistances(0)
        else: 
            for tree in self.trees:
                tree.depthfirst_resistances(0)
            
    def update_forest_tree_volumes(self,indexes, subset):
        if subset:
            for i in indexes:
                self.trees[i].update_volume()
        else:
            for tree in self.trees:
                tree.update_volume()
        
    def test_forest_connection(self, tree_ind, node_index, new_location):       
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
            useless=[]
            self.update_forest_tree_resistances(useless,False) 
            self.update_forest_tree_volumes(useless,False)
            self.node_nber = self.node_nber +2
            return True
        else:
            return False
            
    def printing(self):
        for tree in self.trees:
            print "tree index", tree.tree_index
            tree.printing_full()