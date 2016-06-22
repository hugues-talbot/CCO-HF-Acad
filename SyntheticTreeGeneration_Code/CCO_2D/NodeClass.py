import numpy as np
import copy
import Kamiya as kami

import CCO_2DFunctions as cco_2df

# We define here the classes for our dichotomic tree structure:

# a tree is a list of nodes 
# each node represents a segment extremity, 
# each node has its proper location and index, and contains its parent and two children indexes
# each node also contains information relative to the segment defined with its parent node : flow, resistance, label,beta values
# a tree can be added node one at a time, and has functions to update nodes informations


# this structure is defined this way in order to 
# - avoid information duplicates (that would happen for a tree segment-based)
# - keep manipulations simple: during tree growth we add nodes and modify several nodes information (resistance, beta values, flow),
#   so to run several tests from the same starting tree we decided to produce several copies of it instead of deleting the changes (less risky)
#   also we didn't change the index of the nodes, which helps us keeping track of the growth step  


class Node:
    "Class defining attributes of each segment during synthetic tree generation"
    def __init__(self, index, coord, flow, parent_index):
        self.index = index
        self.coord = coord
        self.flow  = flow
        self.parent_index = parent_index
        self.betas = np.array([1.,1.]) #first if for child_0, second is for child_1
        self.children_index = np.array([-1,-1]) #child_0 index, child_1 index
        self.resistance = -1.
        self.label = -1 # to help debugging
    	
    def set_coord(self, coord):
        self.coord = coord
    		
    def set_flow(self, flow):
        self.flow = flow
    	
    def set_parent_index(self, parent_index):
        self.parent_index = parent_index
    	
    def set_betas(self, betas):
        self.betas = betas
    	
    def set_child_0_index(self, child_0_index):
        self.children_index[0] = child_0_index
    	
    def set_child_1_index(self, child_1_index):
        self.children_index[1] = child_1_index

    def set_resistance(self, resistance):
        print "setting resistance of", self.index, resistance
        self.resistance = resistance
        
    def set_label(self, label):
        self.label = label
        
    def index(self):
        return self.index
    
    def coord(self):
        return self.coord
        
    def betas(self):
        return self.betas
        
    def parent(self):
        return self.parent_index
        
    def children(self):
        return self.children_index
    
    def is_leaf(self):
        if (np.array_equal(self.children_index, np.array([-1,-1]))):
            return True
        else:
            return False
    
    def printing(self):
        print "index", self.index, "parent", self.parent_index, "children", self.children_index, "label", self.label
        print "coord", self.coord, "flow", self.flow, "betas", self.betas, "resistance", self.resistance
        




class Tree:
    "Class defining attributes of the tree during growth"
    def __init__(self, nodes, n_term, q_perf, p_drop, visc):
        self.nodes = nodes
        self.n_term = n_term
        self.q_perf = q_perf
        self.p_drop = p_drop
        self.node_index = 0
        self.nu = visc
    
    def __deepcopy__(self, tree):
        return Tree(copy.deepcopy(self.nodes), copy.deepcopy(self.n_term), copy.deepcopy(self.q_perf), copy.deepcopy(self.p_drop), copy.deepcopy(self.nu))    		
            
    def nodes(self):
        return self.nodes
        
    def update_node_index(self):
        self.node_index = len(self.nodes)
        
    def update_flow(self):
        q_term = self.get_q_term()
        for i in self.nodes:
            if i.is_leaf():
                i.set_flow(q_term)
            else:
                terms = self.get_terms(i.index)
                i.set_flow(terms * q_term)
    
    # get the number of terminal segments in the current tree            
    def get_k_term(self):
        self.k_term = len(self.nodes)/2
        return self.k_term
    		
    def get_q_term(self):
        return float(self.q_perf) / self.get_k_term()   

    def add_node(self, node):		
        if len(self.nodes) < 2*self.n_term:
            self.nodes.append(node)
            self.node_index = self.node_index + 1
            print "node added, total of nodes", self.node_index
            return True
        else:
            print "can't add node, Nterm reached"	
            return False
    		
    def get_node(self, i):
        if (i < len(self.nodes)):		
            return self.nodes[i]
        else :
            print "index exceeds list size"
       
    def printing(self):
        for i in self.nodes:
            print "index", i.index, "coord", i.coord, "betas", i.betas     
        pass
    
    def printing_full(self):
        for i in self.nodes:
            i.printing()
            print " "
        pass
    
    def get_root_radius(self):
        return np.power(self.resistance(1) * self.q_perf / self.p_drop, 1./4)
    
    # get the length of the segment connecting the node of index i and its parent
    def length(self, i):
        parent_ind = self.nodes[i].parent()
        if (parent_ind >= 0):
            p_coord = self.nodes[parent_ind].coord
            i_coord = self.nodes[i].coord
            return np.sqrt(np.sum((p_coord- i_coord)**2))
        else:
            print "no parent found, unable to calculate length"
            return 0.      
              
    def resistance(self, i):
        res = 8*self.nu / np.pi * self.length(i) 
        node = self.get_node(i)
        if i > 0 and (node.is_leaf() == False):
            children_index = node.children()
            res = res + 1./( ( (node.betas[0])**4 / (self.nodes[children_index[0]]).resistance) + 
                          ( (node.betas[1])**4 / (self.nodes[children_index[1]]).resistance) )
        print "total resistance of index ", i, "is ", res
        return res
        
    def update_resistance(self, index):
        print "update resistance of", index
        res = self.resistance(index)
        node = self.nodes[index]
        node.set_resistance(res)
        
    def get_radius(self, index):
        root_radius = self.get_root_radius()
        beta_prod = 1
        while (index > 0):
            parent_ind = (self.nodes[index]).parent()
            parent = self.nodes[parent_ind]            
            if (parent.children()[0] == index):
                beta_prod = beta_prod * parent.betas[0]
            else:
                beta_prod = beta_prod * parent.betas[1]
            index =  parent.index
                    
        return root_radius * beta_prod
          
    
    def volume_iter(self, root_radius, beta_prod, index, volume):
        l = self.length(index)
        volume = volume + l*(beta_prod*root_radius)**2     
        if (self.nodes[index].children())[0] != -1:
            new_beta_prod1 = beta_prod*(self.nodes[index]).betas[0]
            volume = self.volume_iter(root_radius, new_beta_prod1, (self.nodes[index]).children()[0], volume)
        if (self.nodes[index].children())[1] != -1:
            new_beta_prod2 = beta_prod*(self.nodes[index]).betas[1]
            volume = self.volume_iter(root_radius,new_beta_prod2, (self.nodes[index]).children()[1], volume)
        return volume
    
    #light method recursive going following children nodes
    def volume2(self):
        volume_found = 0.        
        root = self.get_root_radius()
        beta_start = 1.
        vol_final = self.volume_iter(root, beta_start, 0, volume_found)
        return vol_final
     
    # beta has been calculated with child_0/child_1 radius ratio where child_0 is old_child         
    def make_connection(self, old_child_index, new_child_location, branching_location, f, beta):
        #creating new nodes
        self.update_node_index()
        parent_index = (self.nodes[old_child_index]).parent()       
        branching_node = Node(self.node_index, branching_location, f[0], parent_index) #flow need to be recalculate everywhere?
        branching_node.set_child_0_index(old_child_index)
        branching_node.set_child_1_index(self.node_index + 1)
        branching_node.set_betas(beta)
        #need to update branching node resistance !!!
        if (self.add_node(branching_node) == False):         
            return False   
            
        new_child_node = Node(self.node_index, new_child_location, f[2], branching_node.index)
        
        length_new_child = np.sqrt(np.sum((branching_location - new_child_location)**2))
        print "length new child", length_new_child
        new_child_resistance = 8.* self.nu * length_new_child / np.pi
        new_child_node.set_resistance(new_child_resistance)
        
        if (self.add_node(new_child_node) == False):
            return False
                    
        #updating the already existing ones:     
        #old child
        old_child_node = self.get_node(old_child_index)
        old_child_node.set_parent_index(branching_node.index)
        old_child_node.set_flow(f[1])
        old_child_node.set_resistance(self.resistance(old_child_index))       
        #branching node 
        branching_node_in_tree = self.nodes[self.node_index-2]
        branching_node_in_tree.set_resistance(self.resistance(self.node_index - 2))       
        #parent
        parent_node = self.get_node(parent_index)
        children_index = parent_node.children()  
        children_index[np.where(children_index == old_child_index)[0]] = branching_node.index
        #no need to update parent_node resistance, will be done in balancing_ratio
        
        return True
                
    def get_terms_recurs(self, index, nb):
        node = self.nodes[index]
        if (node.is_leaf()):
            nb = nb + 1
        else:    
            children_index = node.children()
            if children_index[0] >= 0 :
                nb = self.get_terms_recurs(children_index[0], nb)
            if children_index[1] >= 0 :
                nb = self.get_terms_recurs(children_index[1], nb)
        return nb
        
    # get the number of terminal segments downstream the one of provided index
    def get_terms(self, index):
        nber_of_terms = 0
        nber_of_terms = self.get_terms_recurs(index, nber_of_terms)
        return nber_of_terms
        
    def find_neighbors(self, location, list_size):
        distances=[]
        dist_type = [("distance", float), ("index", int)]
        for i in self.nodes:
            if (i.index > 0):
                dist = cco_2df.segment_distance(i.coord, (self.nodes[i.parent()]).coord, location)
                distances.append((dist, i.index))
        threshold = list_size if (len(self.nodes) > list_size) else len(distances)
        d_array = np.array(distances, dtype = dist_type)
        d_sorted = np.sort(d_array, order = "distance")
        return [i[1] for i in d_sorted[0 : threshold]]  
        
    # test that none of the 3 segments composing the new bifurcation intersect with the other tree segments
    def check_intersection(self, old_child_index, new_child_location, branching_location, new_branches_radii):
        old_child = self.nodes[old_child_index]
        for i in self.nodes:
            if (i.index != old_child_index):
                parent_i = self.nodes[i.parent()]
                radius_i = self.get_radius(i.index)
                print "testing connection with segment", "parent", parent_i.coord, "child", i.coord
                if (cco_2df.no_overlap(i.coord, parent_i.coord, new_child_location, branching_location, radius_i, new_branches_radii[2]) == False):
                    return False
                    
                old_parent_index = old_child.parent()
                old_parent = self.nodes[old_parent_index]
                siblings = old_parent.children()
                old_child_sibling = siblings[0] if (old_child_index == siblings[1]) else siblings[1]
                if (i.index != old_parent_index) and (i.index != old_child_sibling):
                    if (cco_2df.no_overlap(i.coord, parent_i.coord, branching_location, old_parent.coord, radius_i, new_branches_radii[0]) ==  False):
                        return False

                if (old_child.is_leaf() == False):
                    old_child_children = old_child.children()
                    if (i.index != old_child_children[0]) and (i.index != old_child_children[1]):
                        if (cco_2df.no_overlap(i.coord, parent_i.coord, old_child.coord, branching_location, radius_i, new_branches_radii[1]) ==  False):
                            return False
        return True
                   
    def calculate_betas_of_parent(self, index): 
        current_node = self.nodes[index]        
        parent = self.nodes[current_node.parent()]
        print "calculate betas of ", parent.index
        
        if parent.index > 0:
            current_node_is_child_0 = (parent.children()[0] == index)
            sibling_index = parent.children()[1] if current_node_is_child_0 else parent.children()[0]
            sibling = self.nodes[sibling_index]
            resistance_current_child = self.nodes[index].resistance
            resistance_sibling = self.nodes[sibling_index].resistance
            if (current_node_is_child_0):
                sibling_ratio = np.power(current_node.flow * resistance_current_child / (sibling.flow * resistance_sibling), 1./4)                            
            else :
                sibling_ratio = np.power(sibling.flow * resistance_sibling / (current_node.flow * resistance_current_child), 1./4)                      
            betas = cco_2df.calculate_betas(sibling_ratio, 3.)        
            return betas
        else:
            print "beta of root point, no need to calculate"
            return np.array([1., 1.])            

            
    def balancing_ratios(self, index):
        parent_index = self.nodes[index].parent()
        print "balancing ratio of parent of index " , parent_index
        if (parent_index > 0):
            parent = self.nodes[parent_index]            
            betas = self.calculate_betas_of_parent(index)
            parent.set_betas(betas)
            self.update_resistance(parent.index)
            self.balancing_ratios(parent_index)
        else : 
            print "root reached"
                

    # testing the connection between the new_child_location and the segment made of "old_child_index" and its parent   

    
    def test_connection(self, old_child_index, new_child_location):
        # update flow values in tree
        q_term = self.q_perf / (self.get_k_term() + 1) 
        f= np.zeros(3)
        f[1] = self.get_terms(old_child_index) * q_term
        f[2] = q_term
        f[0] = f[1] + f[2]
        
        #set original radii: use the original old child radius for all      
        radius_ori = self.get_radius(old_child_index)
        r = np.ones(3) * radius_ori
        
        #call kamiya: local optimization of the single bifurcation
        iter_max = 100
        tolerance = 0.01
        c0 = (self.nodes[self.nodes[old_child_index].parent()]).coord
        c1 = (self.nodes[old_child_index]).coord
        c2 = new_child_location
        convergence, res = kami.kamiya_single_bif_opt(r, f, c0, c1, c2, iter_max, tolerance, False)
        
        # constraints of the local optimization
        #check boundaries: segments are not degenerated
        if (convergence == True):                     
            branching_location = np.array([res[0][0], res[0][1]]) 
            radii = res[0][2] 
            print "radii", radii
            convergence = cco_2df.degenerating_test(c0,c1,c2, branching_location, radii)
            #check intersection: no segment intersects with other existing tree segments
            if (convergence == True):
                convergence = self.check_intersection(old_child_index, c2, branching_location, radii)
                
        result = [[0.,0.], [0.,0.]]
        if (convergence == True):
            #make connection                       
            estimate_betas = cco_2df.calculate_betas(radii[1]/radii[2], 3.) #this involves old_child will be child_0, new_child is child_1
            
            self.make_connection(old_child_index, new_child_location, branching_location, f, estimate_betas)
            
            #update flows 
            self.update_flow()
                
            # balancing ratio up to the root: update whole tree resistance by tuning radii
            self.balancing_ratios(old_child_index)
            correct_beta= self.nodes[self.node_index-2].betas
            
            #measure tree volume
            tree_volume = self.volume2()
            result=[correct_beta, branching_location]              
            return convergence, tree_volume, result, old_child_index
        else:
            print "connection test failed"
            return False, 0., result, old_child_index
    
    # add the two nodes and update tree
    def add_connection(self, old_child_index, new_child_location, branching_location, betas):
        f= np.zeros(3)
        if (self.make_connection(old_child_index, new_child_location, branching_location, f, betas)):
            self.update_flow()
            self.balancing_ratios(self.node_index - 2)
            print "connection added"
            return True
        else:
            print "no connection added"
            return False
            



            
 

    
