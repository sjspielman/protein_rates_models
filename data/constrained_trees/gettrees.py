import os
from random import uniform
from ete3 import Tree
from copy import deepcopy
import numpy as np

def perturb_bl(etree, fix = None, shape = 0.75):
    if fix is None:
        
        for node in etree.traverse("postorder"):
            if not node.is_root():
                node.dist = np.random.gamma(shape, scale = 1./shape)
    else:
        for node in etree.traverse("postorder"):
            if not node.is_root():
                node.dist = float(fix)
    treestring = etree.write(format=1).strip()

    return treestring
                
                
                
                
                

types = ["enzyme", "virus", "organelle"]

for type in types:
    dir = "../" + type + "_data/"
    files = [x for x in os.listdir(dir) if x.endswith("fna")]
    for file in files:
        
        name = file.split(".fna")[0]
        print name
        with open(dir + file, "r") as f:
            stuff = f.read()
        stuff = stuff.strip().split("\n")
        tree = stuff.pop(-1)

        assert(tree.startswith("(") and tree.endswith(";")), "not a tree"
        etree = Tree(tree, format = 1)    
        
        perttrees = {}
        
        perttrees["fix1"]  = perturb_bl(deepcopy(etree), fix = 0.1)
        perttrees["fix2"]  = perturb_bl(deepcopy(etree), fix = 0.5) ## rates for this tree should be 20% of tree with bl=0.1 if everything is working
        perttrees["runif"] = perturb_bl(deepcopy(etree))

        for key in perttrees:
            with open(type + "_data/" + name + "_" + key + ".fna", "w") as f:
                f.write("\n".join(stuff))
                f.write(perttrees[key])
        assert 1==3

   
