from Bio import AlignIO
from tqdm import tqdm
import numpy as np
import sys
import scipy

try:
    aln = AlignIO.read(sys.argv[1], "fasta")
except:
    aln = AlignIO.read(sys.argv[1], "stockholm")

#AlignIO.write(aln, f"{sys.argv[2]}.aln", "fasta")
dct_q = {}
alphabet = {}
for rec in tqdm(aln, "loading alphabet table"):
    for c in set(rec.seq):
        c = c.lower()
        if c not in alphabet:
            alphabet[c] = len(alphabet)+1

abase = len(alphabet) + 2
nbase = len(alphabet) + 3
gap = -1
for k in alphabet:
    if k == "-":
        gap = alphabet[k]
    elif k == "n":
        alphabet[k] = len(alphabet) + 3
    elif k not in ["a", "t", "c", "g", "-", "n"]:
        alphabet[k] = len(alphabet) + 2

for rec in tqdm(aln, "tranforming seq to vec"):
    arr = np.zeros(shape=(aln.get_alignment_length(),))
    for i, s in enumerate(rec.seq):
        arr[i] = alphabet[s.lower()]
    dct_q[rec.name] = arr

key_list = list(dct_q.keys())
sim_mat = np.zeros(shape=(len(key_list),len(key_list)))
hits = {}
for i in tqdm(range(len(key_list)), "calculating pairwise distance"):
    key = key_list[i]
    qs = dct_q[key]
    hits[key] = []
    for j in range(i+1, len(key_list)):
        d = key_list[j]
        dbs = dct_q[d]
        ide_mask = qs == dbs
        amb_mask = (qs != nbase) & (qs != abase) & (dbs != nbase) & (dbs != abase)
        gap_mask = (qs != gap) & (dbs != gap)
        ident = np.sum(ide_mask & amb_mask & gap_mask)
        aln_len = np.sum(gap_mask & amb_mask)
        pident = ident / aln_len if aln_len > 0 else 0
        hits[key].append((key, d, pident, aln_len))
        sim_mat[i,j] = pident

with open(f"{sys.argv[2]}.tsv", "w") as f:
    for k in hits:
        hits[k].sort(key=lambda x:x[2], reverse=True)
        for query, hit, pident, aln_len in hits[k]:
            f.write(f"{query}\t{hit}\t{pident}\t{aln_len}\n")

