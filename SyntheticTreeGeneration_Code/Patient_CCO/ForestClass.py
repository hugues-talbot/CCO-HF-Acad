# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 14:48:47 2016

@author: jaquetc
"""

import numpy as np
import CCO_3DFunctions as cco_3df

import NodeClass3D as nclass
from scipy.interpolate import RegularGridInterpolator

EPS = 0.001*5 #parameter of tolerance related to the target w value to reach for Newton Algo
MAX_ITER_NEWTON = 50 #max number of iteration for Newton algo
INITIAL_FAC = 0.5
DIST_TO_CENTERLINE = 1.5
LV_INNER_WALL = 2.0
LV_OUTER_WALL = 1.0

class Forest:
    def __init__(self, trees, n_term, q_perf, p_term, visc, v_perf, r_f, vox_size, heart_w, left_a,chamber_w_z, vs_dist, max_curv_radius, im_size, gamma):
        self.trees = trees
        self.n_term = n_term
        self.final_q_perf = q_perf
        self.q_term = q_perf / float(n_term)
        self.p_term = p_term
        self.node_nber = 0 #total of node in forest
        self.nu = visc
        self.v_perf = v_perf
        self.r_supp = np.power(3* v_perf / (4* n_term * np.pi), 1./3) #microbox radius: average radius of terminal segment perfusion territory when reaching final tree growth
        self.final_perf_radius = r_f #real size of the estimated plain perfusion territory radius (excluding concavity volume)
        self.length_factor = 1 # not a const, is updated during tree growth after each added bifurcation       
        self.max_curv_rad = max_curv_radius
        self.im_size = im_size #real size of the final cercle (includes the concavity)
        self.gamma = gamma
        self.voxel_size = vox_size
        self.heart_w = heart_w
        self.empty_grid=(np.arange(0,im_size[0]*vox_size[0],vox_size[0]),np.arange(0,im_size[1]*vox_size[1],vox_size[1]), np.arange(0,im_size[2]*vox_size[2],vox_size[2]))
        self.l_a =  RegularGridInterpolator(self.empty_grid,left_a*4., method ="nearest", bounds_error=False, fill_value = -2.)        
        self.interp_h =   RegularGridInterpolator(self.empty_grid,heart_w + chamber_w_z, bounds_error=False, fill_value = -2.)
        self.gdt_h = np.gradient(heart_w + chamber_w_z)
        self.interp_ghx = RegularGridInterpolator(self.empty_grid,self.gdt_h[0])
        self.interp_ghy = RegularGridInterpolator(self.empty_grid,self.gdt_h[1])
        self.interp_ghz = RegularGridInterpolator(self.empty_grid,self.gdt_h[2])
        self.v_dist_map = RegularGridInterpolator(self.empty_grid,vs_dist, bounds_error=False, fill_value = -2.)
        self.dm_gdt = np.gradient(vs_dist)
        self.dm_gx = RegularGridInterpolator(self.empty_grid, self.dm_gdt[0])
        self.dm_gy = RegularGridInterpolator(self.empty_grid, self.dm_gdt[1])
        self.dm_gz = RegularGridInterpolator(self.empty_grid, self.dm_gdt[2])
        self.dist_max=np.zeros((1,1))
        self.interp_list=[self.interp_h, self.interp_ghx, self.interp_ghy, self.interp_ghz, self.v_dist_map, self.dm_gx, self.dm_gy, self.dm_gz, self.l_a]
    ##trees should be created according to increasing flow 
    def create_tree(self, source_location, vect_dir, q_perf, p_perf, diam_ori):
        tree_index = len(self.trees)
        delta_p = p_perf - self.p_term
        print "delta p of tree", delta_p
        print "tree of index", tree_index
        if q_perf < 20.: #source along vessel
            #change source location
            #if source outside lv use gdt potential vector direction
            #else (inside lv) use inverse direction
            #use diameter/2 distance
            gdt = np.array([self.trees[0].get_hgx(source_location), self.trees[0].get_hgy(source_location), self.trees[0].get_hgz(source_location)])
            gdt_norm_ortho = gdt / cco_3df.length(gdt, self.voxel_size)*self.voxel_size
            pot = self.trees[0].get_h(source_location)
            source_ortho_world = source_location*self.voxel_size
            if pot < LV_OUTER_WALL:
                new_source_ortho = source_ortho_world + 0.5*diam_ori*gdt_norm_ortho
                print "dist between ori and new source", cco_3df.length(source_location-new_source_ortho/self.voxel_size, self.voxel_size) 
                print "0.5*diam", 0.5*diam_ori
                source_location = new_source_ortho / self.voxel_size
            else:               
                new_source_ortho = source_ortho_world - 0.5*diam_ori*gdt_norm_ortho
                print "dist between ori and new source", cco_3df.length(source_location-new_source_ortho/self.voxel_size, self.voxel_size) 
                print "0.5*diam", 0.5*diam_ori
                source_location = new_source_ortho / self.voxel_size
                
        tree = nclass.Tree(tree_index, [], vect_dir, self.q_term, q_perf, delta_p, self.nu, self.v_perf, self.final_perf_radius, self.interp_list, self.voxel_size, self.im_size,self.max_curv_rad, np.max(self.im_size), self.gamma)        
        root_node = nclass.Node(0,source_location, q_perf, -1)
        root_node.set_child_0_index(1)
        tree.add_node(root_node)  
        self.add_tree(tree)     
        
    def add_tree(self, tree):
        self.trees.append(tree)
        self.node_nber = self.node_nber + len(tree.nodes)
        
    def first_segment_end(self, surface_tol):
        dist_max = self.find_sources_dist_max() 
        self.dist_max=dist_max
        #if different flows need to add first segment on biggest flow tree 
        #first forest segment location is constrained to a maximum distance induced by other existing sources location and flow relationships
        indexes = []
        first_added = 0       
        indexes.append(first_added)
        self.trees[first_added].set_activity(True)
#        #testing vectors     
#        for tree in self.trees:
#            node_pos = tree.nodes[0].coord + 10* (tree.direct_vect/ cco_3df.length(tree.direct_vect, self.voxel_size))
#
#            first_node = nclass.Node(1, node_pos, tree.final_q_perf,0)
#            tree.add_node(first_node)
        print "look for first segmt node position"
        first_node_position = self.first_source_end(self.trees[first_added], dist_max[first_added][0]*dist_max[first_added][1],  surface_tol)
        print "dirst node position found!"
        first_node = nclass.Node(1, first_node_position, (self.trees[0]).final_q_perf,0)
        (self.trees[first_added]).add_node(first_node)
        #updates
        (self.trees[first_added]).update_flow()
        self.update_forest_length_factor()
        (self.trees[first_added]).depthfirst_resistances(0)
        (self.trees[first_added]).update_volume()
        self.node_nber = self.node_nber +1        
                            
        for tree in self.trees:
            current_index = tree.tree_index            
            if current_index != first_added:  
                location_found = False
                ind = 0
                while location_found == False:          
                    result = self.get_new_loc_for_sources( 1., tree, dist_max[current_index][0]*dist_max[current_index][1], surface_tol) #1 is dtreshfactor
                    print "tree", tree.tree_index, "h", tree.get_h(result[1])
                    print "dist max", dist_max[current_index], "seg length", cco_3df.length(result[1]-tree.nodes[0].coord, self.voxel_size)
                    location_found = result[0]
                    ind = ind +1
                    print "ind in while loop", ind
                    if ind > 1:
                        print "fail finding location"
                        break
                secund_node = nclass.Node(1, result[1], (self.trees[current_index]).final_q_perf, 0)               
                (self.trees[current_index]).add_node(secund_node)
                #updates
                (self.trees[current_index]).update_flow()
                self.update_forest_length_factor()
                self.update_forest_tree_resistances(indexes, True)
                self.update_forest_tree_volumes(indexes, True)
                self.node_nber = self.node_nber +1
                indexes.append(current_index)  
                print "first seg found for tree of index ", tree.tree_index, result[1]
        self.update_activation()
        #print "dist_max", dist_max      
        print "sources initialized"
                
    def update_activation(self):
        #biggest tree flow has index 0
        biggest_flow = self.trees[0].final_q_perf
        current_flow = self.trees[0].get_q_perf_k()
        for tree in self.trees:
            ratio = biggest_flow / tree.final_q_perf
            if tree.get_q_perf_k() * ratio <= current_flow :
                tree.set_activity(True)
            else:
                tree.set_activity(False)
		 
    
    def new_pos_inside_sqr(self,source_loc, dist_max):
            inf_dist_max = False
            new_pos = np.zeros(2)
            while (inf_dist_max == False):
                rdm_in_sqre = np.random.rand(3)*2.*dist_max
                #print "rdm_in_sqre",rdm_in_sqre
                new_pos = source_loc - (np.ones(3)*dist_max)/self.voxel_size + rdm_in_sqre/self.voxel_size
                print "length in square", cco_3df.length(new_pos-source_loc, self.voxel_size), "dist max", dist_max
                if cco_3df.length(new_pos-source_loc, self.voxel_size) < dist_max:  
                    inf_dist_max = True
            return new_pos
        
    def first_source_end(self,tree, max_dist, surface_tol):
        # the location of the first segment end is constrained by perfusion territory and other sources locations
        source_loc = tree.nodes[0].coord        
        inside_area = False
        while inside_area == False :    
            position = self.new_pos_inside_sqr(source_loc, max_dist)
            if tree.test_vessel_direction(0, position, surface_tol) == False:
                print "wrong vessel direction"
                continue
            if (surface_tol > 0.):
                ins, val = tree.inside_heart(source_loc)
                insi, valu = tree.inside_heart(position)
                if (ins):#if source inside heart and out lv
                    print "source in heart surface point"
                    if (valu >= LV_OUTER_WALL and valu < LV_INNER_WALL):
                        position = self.dist_to_lv_via_sampling(source_loc, position, int(cco_3df.length(source_loc-position, self.voxel_size)*2))
                        if position[0] == 0. and position[1] == 0.:
                            continue
                        n = tree.calculate_sampling(self.max_curv_rad, position, source_loc, LV_OUTER_WALL)
                        print "strp"
                        if tree.sample_and_test(position, source_loc, n, val, False, surface_tol) == True:
                            print "first source end found for surface heart", "h val", tree.get_h(position)
                            return position
                        else:
                            print "sample and test return false"
                    else:
                        continue
                else:#if source already in lv
                    print "source in lv surface point" 
                    if (valu < LV_OUTER_WALL and valu >= 0.):                        
                        n = tree.calculate_sampling(self.max_curv_rad, position, source_loc, LV_INNER_WALL)
                        if tree.sample_and_test(position, source_loc, n, -1., False, surface_tol) == True:
                            print "first source end found for surface lv", tree.get_h(position)
                            return position
                    else:
                        continue
            else: 
                print "surface neg"
                if (tree.inside_perf_territory(position) and self.outside_segmented_vessels(position, DIST_TO_CENTERLINE)):
                    inside, val = self.inside_heart(source_loc)
                    n = tree.calculate_sampling(self.max_curv_rad, position, source_loc, LV_INNER_WALL)
                    if tree.sample_and_test(position, source_loc, n, val, False, surface_tol) == True:
                        print "first source end found"
                        return position
                    else:
                        print "loc outside heart"
                        continue

     
    def get_fk_term(self):
        fk_term = 0
        for tree in self.trees:
            fk_term = fk_term + tree.get_k_term()
        return fk_term
        
    def find_sources_dist_max(self):
        #find closest source for each source
        #store the distance and the ratio of flow between the two
        tree_nb = len(self.trees)
        dist_max = np.ones((tree_nb,2))*np.max(self.im_size)*10.
        #print "dist_max_initial", dist_max
        for tree in self.trees:
            for i in range (tree.tree_index+1, tree_nb):
                dist = cco_3df.length(tree.nodes[0].coord - self.trees[i].nodes[0].coord, self.voxel_size)
                if (dist == 0.):
                    print "dist 0", tree.tree_index, " ", tree.nodes[0].coord, " other ", self.trees[i].tree_index, " ", self.trees[i].nodes[0].coord
                if dist_max[tree.tree_index][0] > dist:
                    dist_max[tree.tree_index][0] = dist
                    tree_q_perf = tree.final_q_perf
                    i_q_perf = self.trees[i].final_q_perf
                    dist_max[tree.tree_index][1] = tree_q_perf / (i_q_perf + tree_q_perf) 
                if dist_max[i][0] > dist:
                    dist_max[i][0] = dist
                    tree_q_perf = tree.final_q_perf
                    i_q_perf = self.trees[i].final_q_perf
                    dist_max[i][1] = i_q_perf / (i_q_perf + tree_q_perf)  
        return dist_max
    
    def get_new_loc_for_sources(self, d_tresh_factor, current_tree, dist_max, surface_tol): 
         #d_tresh is used to constrain distance to existing segments while dist_max is used for distance to other sources
        source_coord = current_tree.nodes[0].coord
        k_term = self.get_fk_term()
        d_tresh, r_pk = cco_3df.calculate_d_tresh_3D(self.r_supp, k_term)
        length_factor = r_pk / self.final_perf_radius
        d_tresh_factorized = d_tresh / length_factor * d_tresh_factor#d_tresh is calculated in the current k_world dimensions
        meet_criteria = False
        ind = 0
        max_it = 500
        print "looking for a position respecting dist_max", dist_max, "and d_tresh", d_tresh_factorized
        #should test : dist_max < d_tresh!
        if (dist_max > d_tresh):
            print "impossible to satisfy both dist_max and d_tresh", dist_max, " ",  d_tresh_factorized
            return False, np.zeros(3), d_tresh_factorized
            
        #detect sources with dist_max < distance to lv: generate an end in direction of lv with dist_max length
        ins, val = self.inside_heart(source_coord)
        if ins: #source is in heart
            print "source in heart",val
            #need to check if lv close enough
            dist_to_lv = self.measure_dist_to_lv(source_coord)
            #print "NEW DIST TO LV", dist_to_lv
            print "dist to lv", dist_to_lv, "distmax", dist_max
            if dist_to_lv *2. > dist_max:
                print "source too far from lv"
                #define the end pointas the one on gdt vector line, at dist_max
                #check its h is lower than source point
                new_location = self.short_segmt_end(source_coord, dist_max, True, current_tree.direct_vect)
                ins_s, val_s = self.inside_heart(source_coord)
                ins, val = self.inside_heart(new_location)   
                print "result of short segmt end", new_location, "h", val, "vals", val_s
                if (val >= val_s):
                    print "short segmt end found", new_location, "dist to source", cco_3df.length(new_location-source_coord, self.voxel_size), "dist max", dist_max                       
                    print "short end", current_tree.get_h(new_location)                    
                    return True, new_location, d_tresh_factorized
                else:
                    print "issue val < val_s"
                    return False, new_location, d_tresh_factorized
            else:
                print "close enough to lv"
        else: #source is in lv
            print "source in lv"
            dist_to_lv = self.measure_dist_to_lv_from_inside(source_coord)
            if (dist_to_lv *2. > dist_max):
                new_location = self.short_segmt_end(source_coord, dist_max, False, current_tree.direct_vect) #should follow tree direction instead(risky)?
                print "weird situation but keep going", dist_to_lv, "dist max", dist_max
                return True, new_location, d_tresh_factorized  
            else:
                print "out of source in lv"
                
                
        while (meet_criteria == False and ind < max_it):        
            new_pos = self.new_pos_inside_sqr(source_coord, dist_max)
            if surface_tol > 0:
                if (current_tree.is_on_surface(new_pos, LV_OUTER_WALL)==False):
                    print "n o s"
                    continue
            if dist_max > DIST_TO_CENTERLINE:
                if self.outside_segmented_vessels(new_pos, DIST_TO_CENTERLINE) == False:
                    print "out 0 ", ind
                    ind = ind + 1
                    continue
            if current_tree.test_vessel_direction(0, new_pos, surface_tol) == False:
                print "wrong direction"
                #ind = ind + 1
                continue
            
            if surface_tol > 0.:#staged growth protocol 
                if np.all(new_pos == source_coord):
                    print "exact same location, skipped!"
                    #ind = ind + 1
                    continue
                else: 
                    n = current_tree.calculate_sampling(self.max_curv_rad, new_pos, source_coord, LV_OUTER_WALL)
                    if current_tree.sample_and_test(new_pos, source_coord, n, val, False, surface_tol) == True:
                        return True, new_pos, d_tresh_factorized
                        
            # not Staged Growth:            
            else:  
                print "surface neg"
                if current_tree.inside_perf_territory(new_pos) : #and self.outside_segmented_vessels(new_pos, DIST_TO_CENTERLINE)
                    respect_criteria = True
                    for tree in self.trees:
                        if tree.tree_index != current_tree.tree_index:
                            if tree.test_dist_criteria(new_pos, d_tresh_factorized,True) == False:
                                respect_criteria = False
                                print "dist criteria non respected", current_tree.tree_index
                                break
                    if (respect_criteria == True):                    
                        ins, val = self.inside_heart(source_coord) #if inside heart (so outside lv: will test sample on heart and lv)
                        n = tree.calculate_sampling(self.max_curv_rad, new_pos, source_coord, LV_INNER_WALL)
                        if tree.sample_and_test(new_pos, source_coord, n, val, False, surface_tol) == True:
                            return True, new_pos, d_tresh_factorized
                        else:
                            print "not found for source", source_coord, "dist max", dist_max, "location", new_pos
                else:
                    print "location out, source coord", source_coord, "dist max", dist_max, "location", new_pos
                    current_tree.inside_perf_territory(source_coord)
                    ind = ind + 1
                    continue
                ind = ind + 1
                if ind == max_it:
                    d_tresh_factorized = 0.9*d_tresh_factorized
                    print "using new value for d_tresh: ", d_tresh_factorized
                    ind = 0
                    
        if ind == max_it:
            if val > LV_OUTER_WALL:
                new_location = self.short_segmt_end(source_coord, dist_max, False, current_tree.direct_vect) #should follow tree direction instead(risky)?
                print "stuck", dist_to_lv, "dist max", dist_max
                return True, new_location, d_tresh_factorized
            else:
                new_location = self.short_segmt_end(source_coord, dist_max, True, current_tree.direct_vect)
                print "stucl out lv", dist_to_lv, "dist max", dist_max
                return True, new_location, d_tresh_factorized
           
        return False, np.zeros(3), d_tresh_factorized
    
    def measure_dist_to_lv(self, source):
        n=40
        seg_end = self.short_segmt_end(source, n, True, np.zeros(3))
        print "short seg emnd found"
        closest_to_lv = self.dist_to_lv_via_sampling(source, seg_end,n*2)
        print "closest point", closest_to_lv                 
        return cco_3df.length(closest_to_lv - source, self.voxel_size)
             
        
    def dist_to_lv_via_sampling(self, start_point,end_point, nsub):
        n=nsub
        p1p2_vec = end_point - start_point
        length_vec = cco_3df.length(p1p2_vec, self.voxel_size)
        #n = int(np.ceil(length_vec/0.1 ))
        previous_val = 0.
        print "n", n, "length_vec", length_vec
        for i in range (int(n)):
            loc = start_point + (i / float(n)) * p1p2_vec
            ins, val= self.inside_heart(loc)
            print "dist to lv via sampling", val
            if val < (LV_OUTER_WALL + 5*EPS) and val > (LV_OUTER_WALL - 5*EPS):
                print "found", val
                return (start_point + (float(i)/n)*p1p2_vec)
            else:
                print "previous_val", previous_val, "val", val
                if (val >= previous_val and val < LV_OUTER_WALL):
                    previous_val = val   
                    print "continue", val
                    continue
                else:
                    if val > LV_OUTER_WALL:
                        print "zooming", val
                        return self.dist_to_lv_via_sampling(start_point + (float(i-1)/n)*p1p2_vec, start_point + (float(i)/n)*p1p2_vec, nsub)   
                    else:
                        print "found weird", val
                        new_end = self.short_segmt_end(start_point + (float(i)/n)*p1p2_vec, length_vec, previous_val<1., np.zeros(3))
                        print "end point", end_point, "new end", new_end
                        return self.dist_to_lv_via_sampling(start_point + (float(i)/n)*p1p2_vec, new_end, nsub)   
               
        print "rror", i, "forest", val, "start", start_point
        return np.zeros(3)
        
    def dist_to_lv_via_sampling_inv(self, start_point, end_point, nsub): #startpoint inside lv
        n=nsub
        p1p2_vec = end_point - start_point
        #
        #n = int(np.ceil(length_vec/0.1 )) 
        previous_val = 2.
        for i in range (n):
            loc = start_point + (i / float(n)) * p1p2_vec
            ins, val= self.inside_heart(loc)
            print "dist to lv via sampling inv", val
            if val < (LV_OUTER_WALL + 5*EPS) and val > (LV_OUTER_WALL - 5*EPS):
                print "found", val
                return (start_point + (float(i)/n)*p1p2_vec)
            else:
                if (val <= previous_val and val > LV_OUTER_WALL):
                    previous_val = val   
                    print "continue", val
                    continue
                else:
                    if val < LV_OUTER_WALL:
                        print "zooming", val
                        return self.dist_to_lv_via_sampling_inv(start_point + (float(i-1)/n)*p1p2_vec, start_point + (float(i)/n)*p1p2_vec, nsub)   
                    else:
#                        print "errroor: previous val:", previous_val, "val", val
#                        return np.zeros(3)
                        print "found weird", val
                        length_vec = cco_3df.length(p1p2_vec, self.voxel_size)
                        print "n", n, "lenght", length_vec, "i", i
                        new_end = self.short_segmt_end(start_point + (float(i)/n)*p1p2_vec, length_vec, previous_val<1., np.zeros(3))
                        print "end point", end_point, "new end", new_end
                        return self.dist_to_lv_via_sampling_inv(start_point + (float(i)/n)*p1p2_vec, new_end, nsub)   
        print "rror", "forest", val, "start", start_point
        return np.zeros(3)
     
    #start from a point inside lv, find the surface point               
    def measure_dist_to_lv_from_inside(self, source):
        n=40
        seg_end = self.short_segmt_end(source, n, False, np.zeros(3))
        #print "measure_dist_to_lv_from_inside"
        closest_to_lv = self.dist_to_lv_via_sampling_inv(source, seg_end, n*2)
        #print "closest point", closest_to_lv                 
        return cco_3df.length(closest_to_lv - source, self.voxel_size)
    
    def get_hgx(self, location):
        return float(self.interp_ghx(location*self.voxel_size))
    
    def get_hgy(self, location):
        return float(self.interp_ghy(location*self.voxel_size))  
        
    def get_hgz(self, location):
        return float(self.interp_ghz(location*self.voxel_size))  
        
    def short_segmt_end(self, source_point, max_dist, gdt_dir, tree_vect):
        gdt = np.array([self.get_hgx(source_point), self.get_hgy(source_point), self.get_hgz(source_point)])
        #gdt = np.array([self.trees[0].get_hgx(source_point), self.trees[0].get_hgy(source_point), self.trees[0].get_hgz(source_point)])
        gdt_norm = gdt / np.sqrt(np.sum(gdt**2))#because in non ortho world cco_3df.length(gdt, self.voxel_size)
        #get the gdt
        #print "short segment end"
        if gdt_dir == False:
            gdt_norm = -gdt_norm

        if np.all(tree_vect== np.zeros(3)) == False:
            gdt_norm = (gdt_norm + tree_vect/np.sqrt(np.sum(tree_vect**2))) / 2.##because in non ortho world cco_3df.length(tree_vect, self.voxel_size
            print "gdt norm", gdt_norm

        #print "gdt", gdt, "gdt norm lengt", cco_3df.length(gdt/length_gdt, self.voxel_size)
        seg_end = source_point + gdt_norm * (max_dist/self.voxel_size)#because in non ortho world
        outside =  (self.interp_h(seg_end*self.voxel_size) < 0.)
        i = 1.2
        if outside:       
            while (outside):               
                seg_end = source_point + gdt_norm* (max_dist/self.voxel_size) *(1./i)#because in non ortho world
                outside = (self.interp_h(seg_end*self.voxel_size) < 0.)
                i = i + 0.1 
        #print "out of while", i
        return seg_end
        #normalize it
        
        
    
    def outside_segmented_vessels(self, location, dist_to_centerline):  
        fmm_val = round(self.trees[0].get_fmm(location),3)
        if (fmm_val> dist_to_centerline): 
                return True
        return False
    
    def inside_perf_territory(self, location):
        return self.trees[0].inside_perf_territory(location)
             
        
    def inside_heart(self,location):
        pot_val = round(float(self.interp_h(location*self.voxel_size)),3)
        if (pot_val> 0.) and (pot_val < LV_OUTER_WALL):
            return True, pot_val
        else: 
            return False, pot_val
#       
    def get_new_location(self, d_tresh_factor, surface_tol):   
        k_term = self.get_fk_term()
        d_tresh, r_pk = cco_3df.calculate_d_tresh_3D(self.r_supp, k_term)
        length_factor = r_pk /  self.final_perf_radius
        d_tresh_factorized = d_tresh / length_factor * d_tresh_factor#d_tresh is calculated in the current k_world dimensions
        meet_criteria = False
        ind = 0
        max_it = 100
        while (meet_criteria == False and ind < max_it):
            point = cco_3df.random_location(self.im_size)
            if surface_tol > 0.:
                #get the point on surface
                ins, val = self.inside_heart(point)
                n_l = 40.
                if val > 0.:
                    
                    if val >= LV_OUTER_WALL and val < LV_INNER_WALL: #inside lv
                        seg_end = self.short_segmt_end(point, n_l, False, np.zeros(3))
                        surf_point = self.dist_to_lv_via_sampling_inv(point, seg_end, int(n_l*2))
                    if val < LV_OUTER_WALL: #inside heart
                         print "inside heart but out lv, skipped"
                         continue
#                        seg_end = self.short_segmt_end(point, n_l, True, np.zeros(3))
#                        surf_point = self.dist_to_lv_via_sampling(point, seg_end)
                    if val >= LV_INNER_WALL:
                        print "in concavity so skipped"
                        continue
                    print "surf_point", surf_point
                    if surf_point[0] == 0. and surf_point[1] == 0.:
                        print "surf point zero"
                        return False, np.zeros(3), d_tresh_factorized
#                    if float(self.l_a(surf_point*self.voxel_size)) > 3.5:
#                        print "on left atrium surface"
#                        continue
                    if (self.outside_segmented_vessels(surf_point, DIST_TO_CENTERLINE)):
                        respect_criteria = True
                        for tree in self.trees:
                            if tree.test_dist_criteria(surf_point, d_tresh_factorized, False) == False:
                                respect_criteria = False
                                break
                        if (respect_criteria == True):
                            print "location on surface found"
                            return True, surf_point, d_tresh_factorized
            else:
                
                if (self.inside_perf_territory(point) == True) and (self.outside_segmented_vessels(point, DIST_TO_CENTERLINE)):
                    respect_criteria = True
                    for tree in self.trees:
                        if tree.test_dist_criteria(point, d_tresh_factorized,False) == False:
                            respect_criteria = False
                            break
                    if (respect_criteria == True):
                        print "location found"
                        print point
                        return True, point, d_tresh_factorized
                else:
                    continue
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
            if self.trees[j].activ: #look only for neighbors of activ trees
                distances = distances + (self.trees[j]).find_neighbors(location, n)  
        print "distances", distances
        threshold = n if (self.node_nber - len(self.trees)> n) else len(distances)
        print "node_nber", self.node_nber, "Ncon ", n, "len(distances) ", len(distances)
        print "threshold", threshold
        d_array = np.array(distances, dtype = dist_type)
        d_sorted = np.sort(d_array, order = "distance")
        return [[i[1],i[2]] for i in d_sorted[0 : threshold]] 
        
    #the length factor is a scaling factor for the segment length
    #from the final coordinates world (dimensions of the final perfusion territory)
    #to the current k world (whose perfusion territory size is increased before each added bifurcation)
    def update_forest_length_factor(self):       
        r_pk = np.sqrt((self.get_fk_term() + 1)* self.r_supp**3) #correspond to the plain perfusion territory (excluding concavity volume)
        print "r_pk",r_pk
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
        
    def test_forest_connection(self, tree_ind, node_index, new_location, surface, surface_tol):       
        current_tree = self.trees[tree_ind]
        print "testing neighbor of tree index", tree_ind
        if surface:
            parent_node_index = (current_tree.nodes[node_index]).parent()
            if current_tree.test_vessel_direction(parent_node_index, new_location, surface_tol) == False:
                print "wrong direction"
                return False, 0.,np.zeros(2), np.zeros(3), tree_ind, node_index, 0.
        reslt = current_tree.test_connection(node_index, new_location, surface_tol)
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
                    print "testing intersection with tree number", tree.tree_index
                    add_rad=np.ones(3)*3.5
#                    if surface_tol > 0.:
#                        add_rad = add_rad + current_tree.mid_point_dist_to_lv(beta_and_bif_location[1], np.array([parent_coord, old_child_coord, new_location]))
                    
                    intersect = tree.check_intersection_external_seg(old_child_coord, parent_coord, new_location, beta_and_bif_location[1], new_radii_rescaled, surface_tol)
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
                return True, total_volume, beta_and_bif_location[0],beta_and_bif_location[1], tree_ind, node_index, current_tree_volume 
            else:
                return False, 0.,np.zeros(2), np.zeros(3), tree_ind, node_index, current_tree_volume
        else:
            return False, 0.,np.zeros(2), np.zeros(3), tree_ind, node_index, current_tree_volume
    
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
            self.update_activation()
            return True
        else:
            return False
            
    def printing(self):
        for tree in self.trees:
            print "tree index", tree.tree_index, "activity", tree.activ
            tree.printing_full()
