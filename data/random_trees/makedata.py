### random trees all around
import dendropy
from Bio import AlignIO
import random
import os

names = [x.split(".fna")[0] for x in os.listdir("../enzyme_data/")]
fasta_directory = "../fasta/"
nreps = 10

for name in names:
    print name
    for i in range(1,nreps+1):
        aln = AlignIO.read(fasta_directory + name + ".fasta", "fasta")
        ntaxa = len(aln)

        tree = dendropy.simulate.treesim.birth_death_tree(2, 0.75, num_extant_tips = ntaxa) # whatever params
        for node in tree.postorder_node_iter():
            node.edge_length = random.random() ## 0-1
        tree.is_rooted = False
        tree.update_bipartitions()
        
        with open(name + "_random_tree_" + str(i) + ".fna", "w") as f:
            for x in range(len(aln)):
                f.write(">T" + str(x+1) + "\n" + str(aln[x].seq) + "\n")
            f.write(tree.as_string(schema="newick").split("[&U] ")[1])
    