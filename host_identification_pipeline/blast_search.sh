query=$1
ref=$2
thread=12
evalue=1e-5

# DIAMOND dbs for ICTV VMR and Virus-Host Database
VMRDB=VMR_MSL38v1.dmnd
VHDB=virushostdb.cds.dmnd

# merge query with reference sequences
if [ $2 ]
then
    seqkit seq "$query" "$ref" > merged_query_and_ref.fa
    query=merged_query_and_ref.fa
fi

# replace all non-standard amino-acid with X
seqkit replace \
    -s \
    -p "[^ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]" \
    -r "X" \
    "$query" \
    > query.renamed.fa
rm "$query"
query=query.renamed.fa

# BLAST against VMR
blast_output="query_vs_vmr.txt"
diamond blastp \
    -q "$query" \
    -d "$VMRDB" \
    -o "$blast_output" \
    -p $thread \
    -e $evalue \
    --very-sensitive

# BLAST against Virus-Host DB
blast_output="query_vs_vhdb.txt"
diamond blastp \
    -q "$query" \
    -d "$VHDB" \
    -o "$blast_output" \
    -p $thread \
    -e $evalue \
    --very-sensitive

# make DIAMOND db for query sequences
diamond makedb \
    --in "$query" \
    -d query.dmnd

# BLAST against self (all-vs-all blast)
blast_output="all_vs_all.txt"
diamond blastp \
    -q "$query" \
    -d query.dmnd \
    -o "$blast_output" \
    -p $thread \
    -e $evalue \
    --very-sensitive

# clean up
rm "$query" query.dmnd
