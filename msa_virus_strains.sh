cd raw_data/viral_genomes_mapping

cd split_by_virus_p80
ls *.fa \
    | sed "s/.fa//g" \
    | parallel -j2 "linsi --maxiterate 1000 --thread 6 --nomemsave {}.fa > {}.aln"

cd split_by_virus_p50
ls *.fa \
    | sed "s/.fa//g" \
    | parallel -j2 "linsi --maxiterate 1000 --thread 6 --nomemsave {}.fa > {}.aln"

cd split_by_virus_p30
ls *.fa \
    | sed "s/.fa//g" \
    | parallel -j2 "linsi --maxiterate 1000 --thread 6 --nomemsave {}.fa > {}.aln"