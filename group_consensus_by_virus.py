from Bio import SeqIO
from pathlib import Path
import re

length_map = {}
name_map = {}
with open("raw_data/viral_genomes/mask_stat.txt") as f:
    lines = f.readlines()[1:]

for line in lines:
    line = line.strip()
    if line:
        sp = line.split()
        name_map[sp[1]] = sp[0]
        length_map[sp[0]] = int(sp[2])

virus_dict = {}
path = Path("raw_data/viral_genomes_mapping/map_to_rep_seqs_nobam")
for fasta_file in path.glob("*.fa"):
    lib_id = re.findall(r"^([A-Z0-9]+)\.rep_seqs", fasta_file.stem).pop()
    for rec in SeqIO.parse(fasta_file, "fasta"):
        virus_name = name_map[rec.id]
        if virus_name not in virus_dict:
            virus_dict[virus_name] = []

        mapped = 0
        mapped += str(rec.seq).count("A")
        mapped += str(rec.seq).count("T")
        mapped += str(rec.seq).count("C")
        mapped += str(rec.seq).count("G")

        length = length_map[virus_name]

        pmapped = mapped*1.0 / length

        rec.id = f"{virus_name}_{lib_id}_len{length}_mapped{mapped}"
        rec.description = ""

        if pmapped > 0.8:
            virus_dict[virus_name].append(rec)

for virus_name in virus_dict:
    if len(virus_dict[virus_name]) > 1:
        SeqIO.write(virus_dict[virus_name],
                    f"raw_data/viral_genomes_mapping/split_by_virus/{virus_name}.fa",
                    "fasta")
