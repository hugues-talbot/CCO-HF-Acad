import numpy as np
import copy
import Kamiya as kami
import CCO_2DFunctions as cco_2df
import sys

GLOBAL_INCREMENT = True


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
    "a_perf: total perfusion area of the final tree"
    "r_f: real radius (mm) of the total perfusion area"

    def __init__(self, nodes, n_term, q_perf, p_drop, visc, a_perf, r_f, gamma):
        self.nodes = nodes
        self.n_term = n_term
        self.final_q_perf = q_perf
        self.q_term = q_perf / float(n_term)
        self.p_drop = p_drop
        self.node_index = 0 #index of the next added node
        self.nu = visc
        self.a_perf = a_perf
        self.r_supp = np.sqrt(a_perf / (n_term * np.pi)) #microbox radius: average radius of terminal segment perfusion territory when reaching final tree growth
        self.final_perf_radius = r_f #real size of the perfusion territory radius when plotting result
        self.length_factor = 1 # not a const, is updated during tree growth after each added bifurcation
        self.gamma = gamma
    def __deepcopy__(self, tree):
        return Tree(copy.deepcopy(self.nodes), self.n_term, self.final_q_perf, self.p_drop, self.nu, self.a_perf, self.final_perf_radius, self.gamma)    		
            
    def nodes(self):
        return self.nodes
    
    #counting the nodes in list to update set the node_index for next added segment    
    def update_node_index(self):
        self.node_index = len(self.nodes)
    
    # updating flow by calculating number of terminal segment downstream each node
    # should be optimized with a depth first traversal path?
    def update_flow(self):
        for i in self.nodes:
            if i.is_leaf():
                i.set_flow(self.q_term)
            else:
                terms = self.get_terms(i.index)
                i.set_flow(terms * self.q_term)
    
    # get the number of terminal segments in the current tree            
    def get_k_term(self):
        self.k_term = len(self.nodes)/2 #because it is a dichotomic tree
        return self.k_term    	
    
    #flow at feeding artery, depending on the number of kterm
    def get_q_perf_k(self):
        return float(self.get_k_term()) * self.q_term
    
    #when adding node on tree need to update the node_index
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
    
    # during tree growth the root segment is modified by new connections, consequently its index changes
    def get_root_index(self):
        children = self.nodes[0].children()
        return children[0] if children[0]>0 else children[1]
    
    #root radius depends on the whole tree resistance, the current Q_perf and the p_drop we are aiming at
    def get_root_radius(self):
        return np.power(self.resistance(self.get_root_index()) * (self.get_q_perf_k()) / self.p_drop, 1./4) 

    #the length factor is a scaling factor for the segment length
    #from the final coordinates world (dimensions of the final perfusion territory)
    #to the current k world (whose perfusion territory size is increased before each added bifurcation)
    def update_length_factor(self):
        r_pk = np.sqrt((self.get_k_term() + 1)* self.r_supp**2)
        self.length_factor =  r_pk / self.final_perf_radius

    
    # get the length of the segment connecting the node of index i and its parent
    def length(self, i):
        parent_ind = self.nodes[i].parent()
        if (parent_ind >= 0):
            p_coord = self.nodes[parent_ind].coord
            i_coord = self.nodes[i].coord
            return np.sqrt(np.sum((p_coord- i_coord)**2)) * self.length_factor
        else:
            print "no parent found, unable to calculate length"
            return 0.      
    
    # resistance of a segment depends on segment length and its children's resistance          
    def resistance(self, i):
        res = 8.*self.nu * self.length(i) / np.pi 
        node = self.get_node(i)
        if i > 0 and (node.is_leaf() == False):
            children_index = node.children()
            res = res + 1./( ( (node.betas[0])**4 / (self.nodes[children_index[0]]).resistance) + 
                          ( (node.betas[1])**4 / (self.nodes[children_index[1]]).resistance) )
        return res
    
    def update_resistance(self, index):
        res = self.resistance(index)
        node = self.nodes[index]
        node.set_resistance(res)

    #update resistances on the tree along a post order depth first traversal
    def depthfirst_resistances(self,index):
        children = (self.nodes[index]).children()
        if children[0] > 0:
            self.depthfirst_resistances(children[0])
        if children[1] > 0:
            self.depthfirst_resistances(children[1])
        if index > 0:
            self.update_resistance(index)

    # the radius corresponds to the product of beta values upstream the node (along root path) with the root radius   
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
    
    #recursive method to calculate the total tree volume
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
    
    #calculate the total volume of the tree: go over all nodes to sum their respective volumes
    def volume(self):
        volume_found = 0.        
        root = self.get_root_radius()
        beta_start = 1.
        vol_final = self.volume_iter(root, beta_start, self.get_root_index(), volume_found)
        return vol_final       
    
    ## test if new location is over d_tresh distance from existing tree segments
    def test_dist_criteria(self, location, d_tresh, area):
        if (cco_2df.belongs_to_area(location, area)):    
            for sgmt in self.nodes:
                if (sgmt.parent() >= 0):
                    dist = cco_2df.segment_distance(sgmt.coord, (self.get_node(sgmt.parent())).coord, location)
                    if (dist < d_tresh):
                        return False
        else:
            return False
        return True      
    
    # a new location is a random location constrained by perfusion territory and distance criterion            
    def get_new_location(self, area_descrpt, n_term):   
        k_term = self.get_k_term()
        d_tresh, r_pk = cco_2df.calculate_d_tresh_2D(self.r_supp, k_term)
        length_factor = r_pk /  area_descrpt[1]
        d_tresh_factorized = d_tresh / length_factor #d_tresh is calculated in the current k_world dimensions
        meet_criteria = False
        ind = 0
        while (meet_criteria == False and ind < 1000):
            point = cco_2df.random_location()
            if (self.test_dist_criteria(point, d_tresh_factorized, area_descrpt)):
                print "location found"
                return True, point, d_tresh_factorized
            ind = ind + 1
            if ind == 1000:
                d_tresh_factorized = 0.9*d_tresh_factorized
                print "using new value: ", d_tresh_factorized
                ind = 0
        return False, np.array([0.,0.]), d_tresh_factorized
     
    # add bifurcation on tree structure, and set the two children new resistances 
    # beta has been calculated with child_0/child_1 radius ratio where child_0 is old_child         
    def make_connection(self, old_child_index, new_child_location, branching_location, f, beta):
        #creating new nodes
        self.update_node_index()
        parent_index = (self.nodes[old_child_index]).parent()       
        branching_node = Node(self.node_index, branching_location, f[0], parent_index) #flow need to be recalculate everywhere?
        branching_node.set_child_0_index(old_child_index)
        branching_node.set_child_1_index(self.node_index + 1)
        branching_node.set_betas(beta)
        if (self.add_node(branching_node) == False):         
            return False   
            
        new_child_node = Node(self.node_index, new_child_location, f[2], branching_node.index)      
        length_new_child = np.sqrt(np.sum((branching_location - new_child_location)**2)) *self.length_factor
        new_child_resistance = 8.* self.nu * length_new_child / np.pi
        new_child_node.set_resistance(new_child_resistance)     
        if (self.add_node(new_child_node) == False):
            return False
                    
        #updating the already existing nodes 
        #(and setting the resistance for the newly created):     
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
    
    #recursive method called to calculate number of terminal segments downstream a node           
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
        
    #find the n closest neighbors to the given location  
    def find_neighbors(self, location, n):
        distances=[]
        dist_type = [("distance", float), ("index", int)]
        for i in self.nodes:
            if (i.index > 0):
                dist = cco_2df.segment_distance(i.coord, (self.nodes[i.parent()]).coord, location)
                distances.append((dist, i.index))
        threshold = n if (len(self.nodes) > n) else len(distances)
        d_array = np.array(distances, dtype = dist_type)
        d_sorted = np.sort(d_array, order = "distance")
        return [i[1] for i in d_sorted[0 : threshold]]  
        
    # test that none of the 3 segments composing the new bifurcation intersect with the other existing tree segments
    # this test will avoid the segments that we are trying to connect to 
    # because they will overlap at bifurcations (old_chil, old_parent and old_child_sibling)
    def check_intersection(self, old_child_index, new_child_location, branching_location, new_branches_radii):
        old_child = self.nodes[old_child_index]
        for i in self.nodes:          
            if (i.index != old_child_index) and (i.index > 0): #also avoid source node (not a segment)
                parent_i = self.nodes[i.parent()]
                inv_length_factor = 1./self.length_factor
                radius_i_rescaled = self.get_radius(i.index) * inv_length_factor
                new_branches_radii_rescaled = new_branches_radii * inv_length_factor
                #print "testing connection with segment", "parent", parent_i.coord, "child", i.coord
                if (cco_2df.no_overlap(i.coord, parent_i.coord, new_child_location, branching_location, radius_i_rescaled, new_branches_radii_rescaled[2]) == False):
                    return False
                    
                old_parent_index = old_child.parent()
                old_parent = self.nodes[old_parent_index]
                siblings = old_parent.children()
                old_child_sibling = siblings[0] if (old_child_index == siblings[1]) else siblings[1]
                if (i.index != old_parent_index) and (i.index != old_child_sibling): 
                    if (cco_2df.no_overlap(i.coord, parent_i.coord, branching_location, old_parent.coord, radius_i_rescaled, new_branches_radii_rescaled[0]) ==  False):
                        return False
                        
                if (old_child.is_leaf() == False):
                    old_child_children = old_child.children()
                    if (i.index != old_child_children[0]) and (i.index != old_child_children[1]):
                        if (cco_2df.no_overlap(i.coord, parent_i.coord, old_child.coord, branching_location, radius_i_rescaled, new_branches_radii_rescaled[1]) ==  False):
                            return False
        return True
          
    #for each node, beta values related to each child are calculated 
    #using the ratios between children of flow*resistance
    def calculate_betas_of_parent(self, index): 
        current_node = self.nodes[index]        
        parent = self.nodes[current_node.parent()]        
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
            betas = cco_2df.calculate_betas(sibling_ratio, self.gamma)        
            return betas
        else:
            #print "beta of root point, no need to calculate"
            return np.array([1., 1.])            

    # the new bifurcation added has impacted on the whole tree resistance:
    # requires radius rescaling, which is propagated upstream until reaching root        
    def balancing_ratios(self, index):
        parent_index = self.nodes[index].parent()
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
        ##update perfusion territory: add a microcirculatory black box
        # by updating the length factor 
        self.update_length_factor()  
        # and updating the resistance on whole tree so all radii are rescaled
        self.depthfirst_resistances(0)          

        # update flow values in tree
        f= np.zeros(3)
        f[1] = self.get_terms(old_child_index) * self.q_term
        f[2] = self.q_term
        f[0] = f[1] + f[2]
        
        #set original radii: use the original old child radius for all      
        radius_ori = self.get_radius(old_child_index)
        r = np.ones(3) * radius_ori
        
        #optimizing loop using Kamyia
        iter_max = 100
        tolerance = 0.01
        c0 = (self.nodes[self.nodes[old_child_index].parent()]).coord
        c1 = (self.nodes[old_child_index]).coord
        c2 = new_child_location        
        iterat = 0
        #kamiya uses a convex average as starting point
        #but schreiner starts from the midpoint of the segment we are connecting with
        x,y = (c0 + c1) / 2. 
           
        #calculate original tree volume or take a bigger approximation (easier because no need of adding segment)
        initial_tree_vol = self.volume() * 10.
        while (iterat < iter_max):
            #call Kamiya : local optimization of the single bifurcation
            conv, x_c, y_c, r_c, l = kami.kamiya_loop_r2(x, y, c0, c1, c2, f, r, self.length_factor, self.nu,self.gamma)
            result = [[0.,0.], [0.,0.]]
            if conv ==  False:
                print "kamyia doesnt converge"
                return False, 0., result, old_child_index
            
            #test degenerating segments
            branching_location = np.array([x_c,y_c])

            #create copy of tree and connect new branch given Kamyia's results
            tree_copy = copy.deepcopy(self)
            tree_copy.update_length_factor()
            tree_copy.depthfirst_resistances(0)
            tree_copy.update_flow()
            estimate_betas = cco_2df.calculate_betas(r_c[1]/r_c[2], self.gamma) 
            #we call it "estimate" because it's calculated with segment radii 
            #(not using flow nor resistance, which is the strict method)
            #the "correct" beta comes from balancing_ratio (resistance update)
            tree_copy.make_connection(old_child_index, new_child_location, branching_location, f, estimate_betas)
            tree_copy.update_flow()
            tree_copy.balancing_ratios(old_child_index)           
            correct_beta= tree_copy.nodes[tree_copy.node_index-2].betas
            
            #calculate total tree volume and volume gradient
            tree_vol = tree_copy.volume()
            new_radii = np.array([tree_copy.get_radius(tree_copy.node_index-2), tree_copy.get_radius(old_child_index), tree_copy.get_radius(tree_copy.node_index-1)])
            if cco_2df.degenerating_test(c0,c1,c2, branching_location, new_radii, self.length_factor) == False :
                    print "degenerated segments"
                    return False, 0., [[0.,0.], [0.,0.]], old_child_index 
            vol_gdt = initial_tree_vol - tree_vol
            if np.abs(vol_gdt) < tolerance:
                #if gradients is under tolerance, and iter_max not reached:                
                #test intersection on the tree not connected to the new branch
                if self.check_intersection(old_child_index, c2, branching_location, new_radii) ==  True:
                    result=[correct_beta, branching_location]  
                    print "connection test succeed!"  
                    return True, tree_vol, result, old_child_index
                else:
                    print "intersection test failed"
                    print tree_copy.length_factor, self.length_factor
                    return False, 0., result, old_child_index
            else:
                #provides Kamiya with current position and radii as new starting point
                initial_tree_vol = tree_vol
                x,y,r = x_c,y_c,new_radii
                iterat = iterat + 1
           
        print "connection test failed : iter max reached"
        return False, 0., result, old_child_index        
        
        
    
    # add the two nodes and update tree
    def add_connection(self, old_child_index, new_child_location, branching_location, betas):
        self.update_length_factor()
        self.depthfirst_resistances(0)
        f= np.zeros(3)
        if (self.make_connection(old_child_index, new_child_location, branching_location, f, betas)):
            self.update_flow()
            self.balancing_ratios(self.node_index - 2)
            #print "connection added"
            return True
        else:
            #print "no connection added"
            return False
            



            
 

    
