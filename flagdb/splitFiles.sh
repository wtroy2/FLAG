# Split formatted_gene2refseq.tsv
# Count the total number of lines in the file
total_lines=$(wc -l < formatted_gene2refseq.tsv)

# Specify the number of parts you want
num_parts=5

# Calculate the number of lines per part
lines_per_part=$((total_lines / num_parts))

# Split the file (skip the header for all but the first file)
split -l $lines_per_part --numeric-suffixes=1 --additional-suffix=.tsv formatted_gene2refseq.tsv formatted_gene2refseq

# Split gene2ensembl.tsv
total_lines=$(wc -l < gene2ensembl.tsv)
num_parts=5
lines_per_part=$((total_lines / num_parts))
split -l $lines_per_part --numeric-suffixes=1 --additional-suffix=.tsv gene2ensembl.tsv gene2ensembl

# Split ncbi_genes.csv
total_lines=$(wc -l < ncbi_genes.csv)
num_parts=5
lines_per_part=$((total_lines / num_parts))
split -l $lines_per_part --numeric-suffixes=1 --additional-suffix=.csv ncbi_genes.csv ncbi_genes
