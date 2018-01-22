"""
    Grab site-specific information from alignments:
        1) Number of total amino acids (aka non-gapped)
        2) Number of unique amino acids
    -----------------------------------
        3) Really this is a proxy for knowing about trajectory and associated exch rates.:
            For example, if the change is REALLY UNLIKELY in WAG, then it will have a higher rate in WAG relative to JC69.
            If it is REALLY LIKELY then, then will have a higher rate in JC69 
            If MODERATELY LIKELY should come out to about the same, about the same rate. Note that this and above are sort of interchangable. I'm "sure" about #1
"""

from Bio import AlignIO,SeqIO, Seq
import os
from pyvolve import Genetics
import numpy as np
from copy import deepcopy

g = Genetics()
AA = g.amino_acids
ZERO = 1e-8
def calc_entropy(f):
    h = -1. * np.sum ( f[f > ZERO] * np.log(f[f > ZERO]) )   
    if h == -0.0:
        h = 0.
    return h 
    
    
    
path = "../data/fasta/"
files = [x for x in os.listdir(path) if x.endswith("fasta")]

outfile = "summarized_data/entropy.csv"
header = "dataset,site,entropy\n"
outrows = {}
for file in files:
    
    rows = []
    name = file.split(".fasta")[0]
    print name
    aln = AlignIO.read(path + file, "fasta")
   
    for i in range(len(aln[0])):
        site = str(i+1)
        column = str( aln[:, i] )         
        
        present_aminos = [x for x in column if x in AA]
        column_f = np.zeros(20)
        for x in present_aminos:
            column_f[ AA.index(x) ] += 1
        column_f /= np.sum(column_f)

        rows.append([site, str(calc_entropy(column_f))])
    
    outrows[name] = rows

with open(outfile, "w") as f:
    f.write(header)
    for dataset in outrows:
        for row in outrows[dataset]:
            f.write(dataset + "," + ",".join(row) + "\n")
    


        
        
    
       




