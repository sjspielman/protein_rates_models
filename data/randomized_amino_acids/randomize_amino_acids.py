"""
    recode amino acids randomly.
    do this with like 
    jukes cantor should give more or less the same rates, UNLESS FREQUENCIES END UP MATTERING? this latter bit might be high.
    
    so, what is making these rates?
    either:
    no difference         matrix
    unknown               frequency biases
    unknown               tree topology
    not much difference   branch lengths
    
    YOU SHOULD ALSO RUN PERMUTED BRANCH LENGTHS, see what the above gives you.
"""
import numpy as np
from Bio import SeqIO
import os
aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


names = [x.split(".fna")[0] for x in os.listdir("../enzyme_data/")]
fasta_directory = "../fasta/"
tree_directory  = "../trees/"
twenty = np.arange(20)
nreps = 10

for name in names:
    print name
    records = list(SeqIO.parse(fasta_directory + name + ".fasta", "fasta"))
    for i in range(1, nreps+1):
    
        np.random.shuffle(twenty)
        new = [aminos[x] for x in twenty]
        swap = dict(zip(aminos, new))

        final = {}
   
        for rec in records:
            sequence = str(rec.seq)
            id = str(rec.id)
        
            newseq = ""
            for aa in sequence:
                try:
                    newseq += swap[aa]
                except:
                    newseq += aa
            final[id] = newseq
    
        outaln = name + "_randomized_aa_" + str(i) + ".fna"
        with open(tree_directory + name + ".tre", "r") as f:
            tree = f.read().strip()
        with open(outaln, "w") as f:
            for id in final:
                f.write(">" + id + "\n" + final[id] + "\n")
            f.write(tree)
            
    
    
        
    
   
   

    