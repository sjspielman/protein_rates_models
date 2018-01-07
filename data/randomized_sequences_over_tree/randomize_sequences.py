from random import shuffle
from Bio import SeqIO
import os



names = [x.split(".fna")[0] for x in os.listdir("../enzyme_data/")]
fasta_directory = "../fasta/"
tree_directory  = "../trees/"
nreps = 10

for name in names:
    print name
    aln = list(SeqIO.parse(fasta_directory + name + ".fasta", "fasta"))
    ids = [str(record.id) for record in aln]

    with open("../trees/" + name + ".tre", "r") as f:
        treestring = f.read().strip()
    
    for i in range(1, nreps+1):
    
        outfile = name + "_randomized_sequences_" + str(i) + ".fna"
        shuffle(ids)
    
        with open(outfile, "w") as f:
            for x in range(len(aln)):
                f.write(">" + str(ids[x]) + "\n" + str(aln[x].seq) + "\n")
            f.write(treestring)
    
        