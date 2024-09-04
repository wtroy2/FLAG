cd flagdb

echo "ungzipping files"
gunzip *

echo "putting split files back together"
num_parts=6

# Append the content of each part, skipping headers in parts after the first
for i in $(seq 1 $num_parts); do
    cat gene2ensembl0${i}.tsv >> gene2ensembl.tsv
    cat ncbi_genes0${i}.csv >> ncbi_genes.csv
    rm gene2ensembl0${i}.tsv
    rm ncbi_genes0${i}.csv
done
num_parts=5
for i in $(seq 1 $num_parts); do
    cat formatted_gene2refseq0${i}.tsv >> formatted_gene2refseq.tsv
    rm formatted_gene2refseq0${i}.tsv
done

cd ..

echo "Tarring it up for use for the pipeline"
tar -zcvf flagdb.tar.gz flagdb

echo "Done"