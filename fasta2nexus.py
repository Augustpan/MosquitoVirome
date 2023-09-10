from Bio import SeqIO
import sys

records = list(SeqIO.parse(f"{sys.argv[1]}", "fasta"))
for rec in records:
    rec.annotations["molecule_type"] = "DNA"
count = SeqIO.write(records, f"{sys.argv[2]}", "nexus")
