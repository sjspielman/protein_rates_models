from Bio import AlignIO
import os
import numpy as np


names = ["13PK_A" ,"1C0K_A" ,"1E6E_A" ,"1GP1_A" ,"1MQW_A" ,"1REQ_A" ,"1YON_A" ,"1A16_A" ,"1C2T_A" ,"1E7Q_A"]
aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

for name in names:
    aln = AlignIO.read("../../data/fasta/" + name + ".fasta", "fasta")

    aacounts = np.zeros(20)
    total = 0.
    for record in aln:
        seq = str(record.seq)
        for aa in seq:
            if aa in aminos:
                aacounts[aminos.index(aa)] += 1.
                total += 1.
    aacounts/=total  
    assert((1. - np.sum(aacounts))<=1e-8), "bad f sum"

    new_order = [0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17] 
    # Save to raxml model file
    with open("JC69_" + name + ".raxmat", "w") as raxf:
        for i in range(20):
            for j in range(20):
                if (i==j):
                    rate ="0.0"
                else:
                    rate = "1.0"
                raxf.write(rate + "\n")
        for i in range(20):
            raxf.write(str(aacounts[new_order[i]]) + "\n")
