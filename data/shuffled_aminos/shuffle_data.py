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
aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


name = "../data/fasta/1A16_A.fasta"
records = list(SeqIO.parse(name, "fasta"))
twenty = np.arange(20)

for i in range(10):
    
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
    
    outaln = "1A16_A_shuffledAA_" + str(i) + ".fasta"
    with open(outaln, "w") as f:
        for id in final:
            f.write(">" + id + "\n" + final[id] + "\n")
    
    
    
        
    
   
   

    