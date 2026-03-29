import sys
from Bio import AlignIO
import pandas as pd

aln = AlignIO.read(sys.argv[1], "fasta")

names = [rec.id for rec in aln]
matrix = pd.DataFrame(index=names, columns=names)

for i in range(len(aln)):
    for j in range(len(aln)):
        dist = sum(a != b for a, b in zip(aln[i].seq, aln[j].seq))
        matrix.iloc[i, j] = dist

matrix.to_csv(sys.argv[2])
