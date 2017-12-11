"""
    Translate, tree, align the holmes data.
    Some of the data is bad, so it gets continue'd
"""

from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import subprocess

indir = "oxford/"
outdir = "oxford_processed/"
raw = [x for x in os.listdir(indir) if x.endswith(".nex")]

for file in raw:
    
    name = file.split("_.nex")[0]
    print name
    outf2 = outdir + name + "_AA.aln"
    outf3 = outdir + name + "_AA.tre"   
    outf4 = outdir + name + "_AA.fna"
    if os.path.exists(outf4):
        continue
    
    keep = True
    try:
        records = list( SeqIO.parse(indir + file, "fasta"))
    except:
        continue
    if records == []:
        continue
        
    prot_dict = {}
    for entry in records:
        seq = str(entry.seq).replace("-","")
        if len(seq)%3 != 0 or len(seq) == 0:
            keep=False
            break
        prot_seq = ''
        for i in range(0,len(seq),3):
            codon = seq[i:i+3]
            try:
                amino = str( Seq.Seq(codon, generic_dna).translate() )
            except:
                raise AssertionError("\n\nCould not translate input nucleotide codon, quitting.")
            if amino == '*':
                # If stop codon is the last codon, just remove it. If internal stop codon, freak out!
                if i == len(nucseq)-3:
                    nuc_dict[entry] = nuc_dict[entry][:-3]
                else:
                    raise AssertionError("\n\n You have internal stop codons, quitting.")
            else:
                prot_seq += amino
        prot_dict[str(entry.id)] = prot_seq
    
    if keep is False:
        continue
    outf1 = outdir + name + "_AA.fasta"
    with open(outf1, "w") as f:
        for entry in prot_dict:
            f.write(">" + str(entry) + "\n" + str(prot_dict[entry]) + "\n")

    x = subprocess.call("mafft --auto --quiet " + outf1 + " > " + outf2, shell=True)
    assert(x==0)
    
    x = subprocess.call("FastTree -lg -quiet -nosupport " + outf2 + " > " + outf3, shell=True)
    assert(x==0)    

    x=subprocess.call("cat "  + outf2 + " " + outf3 + " > " + outf4, shell=True)
    assert(x==0)    
    
