from Bio import SeqIO
import sys

records = SeqIO.parse(f"{sys.argv[1]}", "fasta")
count = SeqIO.write(records, f"{sys.argv[2]}", "phylip-relaxed")
