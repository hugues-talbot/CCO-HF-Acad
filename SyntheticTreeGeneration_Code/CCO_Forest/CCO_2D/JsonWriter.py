# -*- coding: utf-8 -*-
"""
Created on Wed Nov 09 10:51:12 2016

@author: cjaquet
"""

#import ForestClass as FC
import json
import numpy as np
def write_json(forest, filename):
    file_out = open(filename, 'w')
    fdict={}
    fdict.update({"Flow" : forest.final_q_perf})
    fdict.update({"Trees" : len(forest.trees)})
    fdict.update({"Nterm" : forest.n_term})
    fdict.update({"FKterm" : forest.get_fk_term()})
    fdict["Trees"]= {}
    for tree in forest.trees:
        fdict["Trees"][tree.tree_index]={}
        fdict["Trees"][tree.tree_index].update({"Tree flow" : tree.final_q_perf})
        fdict["Trees"][tree.tree_index].update({"Tree pressure drop" : tree.p_drop})
        fdict["Trees"][tree.tree_index].update({"Kterm" : tree.get_k_term()})
        fdict["Trees"][tree.tree_index].update({"Resistance": tree.resistance(tree.get_root_index())})
        fdict["Trees"][tree.tree_index]["Nodes"]= {}
        for node in tree.nodes:
            fdict["Trees"][tree.tree_index]["Nodes"][node.index]={}
            coords = node.coord.tolist()
            print "coords",coords
            coords.append(0.)
            print "coords",coords
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Location" : coords})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Parent" : node.parent_index})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Children" : node.children_index.tolist()})
            fdict["Trees"][tree.tree_index]["Nodes"][node.index].update({"Radius" : tree.get_radius(node.index)})
            
    json.dump(fdict, file_out, indent = 4, sort_keys = True)
    file_out.close()
            
    
    return True
    
def write_tree_json(tree, filename):
    file_out = open(filename, 'w')
    fdict={}
    fdict["Trees"]= {}
    tree_index = 0
    fdict["Trees"][tree_index]={}
    fdict["Trees"][tree_index].update({"Tree flow" : tree.final_q_perf})
    fdict["Trees"][tree_index].update({"Tree pressure drop" : tree.p_drop})
    fdict["Trees"][tree_index].update({"Kterm" : tree.get_k_term()})
    fdict["Trees"][tree_index].update({"Resistance": tree.resistance(tree.get_root_index())})
    fdict["Trees"][tree_index]["Nodes"]= {}
    for node in tree.nodes:
        fdict["Trees"][tree_index]["Nodes"][node.index]={}
        coords = node.coord.tolist()
        print "coords",coords
        #coords.append(0.)
        print "coords",coords
        fdict["Trees"][tree_index]["Nodes"][node.index].update({"Location" : coords})
        fdict["Trees"][tree_index]["Nodes"][node.index].update({"Parent" : node.parent_index})
        fdict["Trees"][tree_index]["Nodes"][node.index].update({"Children" : node.children_index.tolist()})
        fdict["Trees"][tree_index]["Nodes"][node.index].update({"Radius" : tree.get_radius(node.index)})
            
    json.dump(fdict, file_out, indent = 4, sort_keys = True)
    file_out.close()
            
    
    return True