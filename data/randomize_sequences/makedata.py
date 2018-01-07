### reassign sequences across tips randomly, 10 replicates

import os
from Bio import SeqIO
from random import shuffle


aln = list(SeqIO.parse("1A16_A.fasta", "fasta")) ### 301
ids = [str(record.id) for record in aln]


with open("../trees/1A16_A.tre", "r") as f:
    treestring = f.read().strip()

for i in range(1,11):
    outfile = "1A16_A_randomizeseqs_" + str(i) + ".fna"
    shuffle(ids)
    
    with open(outfile, "w") as f:
        for x in range(len(aln)):
            f.write(">" + str(ids[x]) + "\n" + str(aln[x].seq) + "\n")
        f.write(treestring)
    
        