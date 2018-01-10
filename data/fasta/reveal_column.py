from Bio import AlignIO
import sys

file = sys.argv[1] + ".fasta"
site = int(sys.argv[2])-1

alignment = AlignIO.read(file, "fasta")

print(alignment[:,site])
