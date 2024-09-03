#!/usr/bin/env python
# coding: utf-8
import csv
from collections import defaultdict
import argparse
from datetime import datetime
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
parser = argparse.ArgumentParser(description="Process gene annotation and related files.")

parser.add_argument('--annotation_file', default=None, help='Path to the annotation GFF3 file')
parser.add_argument('--ncbi_ortholog_file', default='ncbi_orthologs.tsv', help='Path to the NCBI ortholog file')
parser.add_argument('--ensembl_ortholog_file', default='ensembl_orthologs.csv', help='Path to the Ensembl ortholog file')

parser.add_argument('--gene2refseq_file', default='formatted_gene2refseq.tsv', help='Path to the gene2refseq file')
parser.add_argument('--gene2ensembl_mapping_file', default='ensembl_genes.csv', help='Path to the gene2ensembl mapping file')

parser.add_argument('--mapping_file_path', default='gene2ensembl.tsv', help='Path to the gene2ensembl mapping file')
parser.add_argument('--ncbi_ortholog_file_path', default='ncbi_orthologs.tsv', help='Path to the NCBI ortholog file path')
parser.add_argument('--ncbi_exon_mapping_file', default='ncbi_genes.csv', help='Path to the NCBI exon mapping file')

parser.add_argument('--splign_file', default=None, help='Path to the Splign GFF3 file')
parser.add_argument('--reference_splign_file', default=None, help='Path to the reference Splign GFF3 file')
parser.add_argument('--pasa_file', default=None, help='Path to the PASA assembly GFF3 file')
parser.add_argument('--reference_pasa_file', default=None, help='Path to the reference PASA assembly GFF3 file')
parser.add_argument('--busco_file', default=None, help='Path to the BUSCO full table file')

parser.add_argument('--protein_alignment_file', default=None, help='Path to the protein alignment file')

parser.add_argument('--output_evidence_gff', default='annotation_with_evidence.gff3', help='Output GFF3 file with evidence annotations')
parser.add_argument('--output_removed_gff', default='removed.gff3', help='Output GFF3 file for removed annotations')
parser.add_argument('--csv_file_path', default="FLAG_Gene_Support_Report.csv", help='CSV file path for the gene support report')
parser.add_argument('--output_supported_gff', default='supported.gff3', help='Output GFF3 file for supported annotations')

# Arguments for multi exon conditions
parser.add_argument("--multi_percent_same_exon_match", type=int, default=30, help="Minimum percent for same exon match in multi-exon conditions")
parser.add_argument("--multi_splign_coverage_min", type=int, default=75, help="Minimum Splign coverage for multi-exon conditions")
parser.add_argument("--multi_reference_splign_coverage_min", type=int, default=65, help="Minimum reference Splign coverage for multi-exon conditions")
parser.add_argument("--multi_pasa_coverage_min", type=int, default=75, help="Minimum PASA coverage for multi-exon conditions")
parser.add_argument("--multi_reference_pasa_coverage_min", type=int, default=65, help="Minimum reference PASA coverage for multi-exon conditions")
parser.add_argument("--multi_miniprot_coverage_min", type=int, default=85, help="Minimum miniprot coverage for multi-exon conditions")
parser.add_argument("--multi_is_busco", type=bool, default=True, help="Is BUSCO required for multi-exon conditions")
parser.add_argument("--multi_num_orthologs", type=int, default=1, help="Number of orthologs for multi-exon conditions")

# Arguments for single exon conditions
parser.add_argument("--single_percent_same_exon_match", type=int, default=10, help="Minimum percent for same exon match in single-exon conditions")
parser.add_argument("--single_splign_coverage_min", type=int, default=85, help="Minimum Splign coverage for single-exon conditions")
parser.add_argument("--single_reference_splign_coverage_min", type=int, default=75, help="Minimum reference Splign coverage for single-exon conditions")
parser.add_argument("--single_pasa_coverage_min", type=int, default=85, help="Minimum PASA coverage for single-exon conditions")
parser.add_argument("--single_reference_pasa_coverage_min", type=int, default=75, help="Minimum reference PASA coverage for single-exon conditions")
# If miniprot_coverage_min is not applicable for single exon, you can omit or set a default value
parser.add_argument("--single_miniprot_coverage_min", type=int, help="Minimum miniprot coverage for single-exon conditions")
parser.add_argument("--single_is_busco", type=bool, default=True, help="Is BUSCO required for single-exon conditions")
parser.add_argument("--single_num_orthologs", type=int, default=5, help="Number of orthologs for single-exon conditions")
parser.add_argument("--target_ratio_single_multi", type=float, default=0.12, help="Target ratio of the number of single exon genes to multi exon genes. Default is 0.12")
parser.add_argument("--length_threshold", type=int, default=10000000, help="Max gene length")
# multi_exon_conditions = {
#     "percent_same_exon_match": 30,
#     "splign_coverage_min": 75,
#     "reference_splign_coverage_min": 65,
#     "pasa_coverage_min": 75,
#     "reference_pasa_coverage_min": 65,
#     "miniprot_coverage_min": 85,
#     "is_busco": True,
#     "num_orthologs": 1
# }
# single_exon_conditions = {
#     "percent_same_exon_match": 10,
#     "splign_coverage_min": 85,
#     "reference_splign_coverage_min": 75,
#     "pasa_coverage_min": 85,
#     "reference_pasa_coverage_min": 75,
#     # "miniprot_coverage_min": 101,
#     "is_busco": True,
#     "num_orthologs": 5
# }
# Constructing the conditions dictionaries

args = parser.parse_args()

# Assigning each parsed argument to the original variable names
annotation_file = args.annotation_file
ncbi_ortholog_file = args.ncbi_ortholog_file
ensembl_ortholog_file = args.ensembl_ortholog_file

gene2refseq_file = args.gene2refseq_file
gene2ensembl_mapping_file = args.gene2ensembl_mapping_file

mapping_file_path = args.mapping_file_path
ensembl_ortholog_file_path = ensembl_ortholog_file
ncbi_ortholog_file_path = args.ncbi_ortholog_file_path
ncbi_exon_mapping_file = args.ncbi_exon_mapping_file

splign_file = args.splign_file
reference_splign_file = args.reference_splign_file
pasa_file = args.pasa_file
reference_pasa_file = args.reference_pasa_file
busco_file = args.busco_file

protein_alignment_file = args.protein_alignment_file

output_evidence_gff = args.output_evidence_gff
output_removed_gff = args.output_removed_gff
csv_file_path = args.csv_file_path
output_supported_gff = args.output_supported_gff
length_threshold = args.length_threshold

target_ratio_single_multi = args.target_ratio_single_multi
if args.single_miniprot_coverage_min:
    single_exon_conditions = {
        "percent_same_exon_match": args.single_percent_same_exon_match,
        "splign_coverage_min": args.single_splign_coverage_min,
        "reference_splign_coverage_min": args.single_reference_splign_coverage_min,
        "pasa_coverage_min": args.single_pasa_coverage_min,
        "reference_pasa_coverage_min": args.single_reference_pasa_coverage_min,
        "miniprot_coverage_min": args.single_miniprot_coverage_min,
        "is_busco": args.single_is_busco,
        "num_orthologs": args.single_num_orthologs,
    }
else:
    single_exon_conditions = {
        "percent_same_exon_match": args.single_percent_same_exon_match,
        "splign_coverage_min": args.single_splign_coverage_min,
        "reference_splign_coverage_min": args.single_reference_splign_coverage_min,
        "pasa_coverage_min": args.single_pasa_coverage_min,
        "reference_pasa_coverage_min": args.single_reference_pasa_coverage_min,
        "is_busco": args.single_is_busco,
        "num_orthologs": args.single_num_orthologs,
    }

multi_exon_conditions = {
    "percent_same_exon_match": args.multi_percent_same_exon_match,
    "splign_coverage_min": args.multi_splign_coverage_min,
    "reference_splign_coverage_min": args.multi_reference_splign_coverage_min,
    "pasa_coverage_min": args.multi_pasa_coverage_min,
    "reference_pasa_coverage_min": args.multi_reference_pasa_coverage_min,
    "miniprot_coverage_min": args.multi_miniprot_coverage_min,
    "is_busco": args.multi_is_busco,
    "num_orthologs": args.multi_num_orthologs,
}

if annotation_file == None:
    print("This needs an annotation file in gff3 format as input")
    exit(1)


# annotation_file = 'finalTin_gut.gff3'
# ncbi_ortholog_file = 'ncbi_orthologs.tsv'
# ensembl_ortholog_file = 'ensembl_orthologs.csv'

# gene2refseq_file = 'formatted_gene2refseq.tsv'


# gene2ensembl_mapping_file = 'parse_combined_ensembl.csv'

# mapping_file_path = 'gene2ensembl.tsv'

# ensembl_ortholog_file_path = 'ensembl_orthologs.csv'

# ncbi_ortholog_file_path = 'ncbi_orthologs.tsv'

# ncbi_exon_mapping_file = 'parse_combined_filtered_file_ncbi.csv'



# splign_file = 'splign.gff3'
# reference_splign_file = 'reference_splign.gff3'
# pasa_file = 'pasa_assembly.gff3'
# reference_pasa_file = 'reference_pasa_assembly.gff3'
# busco_file = 'full_table.tsv'

# protein_alignment_file = "miniprot.gtf"

# output_evidence_gff = 'annotation_with_evidence.gff3'
# output_removed_gff = 'removed.gff3'
# csv_file_path = "FLAG_Gene_Support_Report.csv"
# output_supported_gff = 'supported.gff3'


def parse_annotation_file(annotation_file):
    """
    Parses the annotation file to extract eggnog_orthologs along with their gene IDs for gene features only.
    Returns a list of unique tuples, each containing the gene ID, eggnog_ortholog, its source (NCBI or Ensembl),
    and taxid (None for Ensembl), preserving their order of appearance.
    """
    seen = set()
    eggnog_orthologs = []
    with open(annotation_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) < 9 or row[2] != "gene":  # Check if the feature type is "gene"
                continue
            attributes = {field.split('=')[0]: field.split('=')[1] for field in row[8].split(';') if '=' in field}
            eggnog_ortholog = attributes.get('eggnog_ortholog', '')
            gene_id = attributes.get('ID', '').split(';')[0]  # Assuming ID is always present and parsing the gene ID
            if eggnog_ortholog != '""':
                if '.' in eggnog_ortholog:
                    parts = eggnog_ortholog.split('.')
                    taxid = parts[0] if parts[0].isdigit() else None
                    eggnog_ortholog = parts[1] if taxid else parts[0]  # Adjusted to ensure correct eggnog_ortholog assignment
                    if 'ENS' in eggnog_ortholog:
                        source = 'Ensembl'
                    else:
                        source = 'NCBI'
                tuple_entry = (gene_id, eggnog_ortholog, source, taxid)
                if tuple_entry not in seen:
                    seen.add(tuple_entry)
                    eggnog_orthologs.append(tuple_entry)
    return eggnog_orthologs


# In[2]:


# Step 2: Parse the gene2ensembl.tsv file
def parse_gene2ensembl(file_path):
    protein_to_gene_ensembl = defaultdict(list)
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            protein_accession_version = row['protein_accession.version'].split('.')[0]  # Splitting to get accession without version
            row_taxid = row['#tax_id']
            # Check if (protein_accession, taxid) tuple is in our set
            if (protein_accession_version, row_taxid) in ncbi_proteins:
                # Add all versions for comparison later, including taxid for verification
                protein_to_gene_ensembl[(protein_accession_version, row_taxid)].append(row)
    return protein_to_gene_ensembl

# Step 3: Select the highest version for each protein accession, ensuring taxid matches
def get_highest_version(protein_to_gene_ensembl):
    highest_version_mapping = {}
    for (protein, taxid), entries in protein_to_gene_ensembl.items():
        highest_version_entry = max(entries, key=lambda x: float(x['protein_accession.version'].split('.')[1]))
        highest_version_mapping[(protein, taxid)] = (highest_version_entry['GeneID'], highest_version_entry['Ensembl_gene_identifier'])
    return highest_version_mapping

def load_orthologs(ortholog_file, source):
    """
    Loads ortholog data from a given file, filtering by source (NCBI or Ensembl).
    Returns a dictionary mapping from gene/protein IDs to sets of orthologous IDs.
    """
    orthologs = {}
    with open(ortholog_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header
        for line in file:
            row = line.strip().split('\t')  # Assuming tab-delimited; change as necessary
            # Skip empty lines or lines with fewer than expected columns
            if not row or len(row) < 3:
                continue
            if source == 'Ensembl':
                # Assuming Ensembl format
                gene_id = row[2]  # Adjust index as per actual column position
            elif source == 'NCBI':
                # Assuming NCBI format
                gene_id = row[1]  # Adjust index as per actual column position
            if source == 'NCBI' and gene_id not in orthologs:
                orthologs[gene_id] = set()
            if source == 'Ensembl' and ensembl_id not in orthologs:
                orthologs[ensembl_id] = set()
            orthologs[gene_id if source == 'NCBI' else ensembl_id].add(row[1] if source == 'Ensembl' else row[2])
    return orthologs


# In[3]:


def load_gene2refseq(gene2refseq_file):
    gene2refseq_mappings = []  # Use a list to store tuples
    with open(gene2refseq_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            tax_id, gene_id, protein_accession = line.strip().split('\t')
            # Only append if protein_accession is not '-'
            if protein_accession != '-':
                gene2refseq_mappings.append((tax_id, protein_accession, gene_id))
    return gene2refseq_mappings


# In[4]:


# Example usage

eggnog_orthologs = parse_annotation_file(annotation_file)
# print(eggnog_orthologs)
# Step 1: Extract relevant NCBI protein accessions and keep track of Form Bio gene IDs for them
ncbi_proteins = {(protein_accession, taxid) for _, protein_accession, source, taxid in eggnog_orthologs if source == "NCBI"}


# In[5]:


gene2refseq_mapping = load_gene2refseq(gene2refseq_file)

print("parsed gene2refseq")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
# # Note not all ensembl proteins have a matching geneid. for this case just say supported by at least 1 ortholog

# In[6]:

def preprocess_gene2ensembl_mapping(filepath):
    lookup_dict = {}
    with open(filepath, 'r') as file:
        next(file)  # Skip header if it exists
        for line in file:
            columns = line.strip().split(',')
            gene_id, protein_id, exon_count, assembly_id = columns
            lookup_dict[protein_id] = gene_id  # Assuming protein_id is unique and directly maps to gene_id
    return lookup_dict


# In[7]:


def preprocess_gene2refseq_mapping(gene2refseq_mapping):
    lookup_dict = {}
    for tax_id, protein_accession, gene_id in gene2refseq_mapping:
        # Remove version from protein_accession
        protein_accession_no_version = protein_accession.rsplit('.', 1)[0]
        lookup_key = (tax_id, protein_accession_no_version)
        lookup_dict[lookup_key] = gene_id
    return lookup_dict

def match_eggnog_with_mappings(eggnog_orthologs, gene2refseq_mapping, gene2ensembl_mapping_file):
    refseq_lookup_dict = preprocess_gene2refseq_mapping(gene2refseq_mapping)
    ensembl_lookup_dict = preprocess_gene2ensembl_mapping(gene2ensembl_mapping_file)
    matched_entries = []

    for formbioid, protein_accession, source_type, taxid in eggnog_orthologs:
        protein_accession_no_version = protein_accession.rsplit('.', 1)[0]
        gene_id = None  # Default to None if no match is found
        if source_type == 'NCBI':
            lookup_key = (taxid, protein_accession_no_version)
            gene_id = refseq_lookup_dict.get(lookup_key, None)  # Use .get to return None if key is not found
            # if lookup_key in refseq_lookup_dict:
            #     gene_id = refseq_lookup_dict[lookup_key]
            #     matched_entries.append((formbioid, protein_accession, source_type, taxid, gene_id))
        elif source_type == 'Ensembl':
            gene_id = ensembl_lookup_dict.get(protein_accession_no_version, None)  # Adjusted to use .get method
            # if protein_accession in ensembl_lookup_dict:
            #     gene_id = ensembl_lookup_dict[protein_accession]
            #     matched_entries.append((formbioid, protein_accession, source_type, taxid, gene_id))
        matched_entries.append((formbioid, protein_accession, source_type, taxid, gene_id))

    return matched_entries


# In[8]:


matched_entries = match_eggnog_with_mappings(eggnog_orthologs, gene2refseq_mapping, gene2ensembl_mapping_file)


# In[9]:


def parse_geneid_mapping(file_path):
    ncbi_to_ensembl = {}
    ensembl_to_ncbi = {}
    with open(file_path, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split('\t')
            tax_id, ncbi_geneid, ensembl_geneid = parts[:3]
            ncbi_to_ensembl[ncbi_geneid] = ensembl_geneid
            ensembl_to_ncbi[ensembl_geneid] = ncbi_geneid
    return ncbi_to_ensembl, ensembl_to_ncbi
def enhance_matched_entries(matched_entries, ncbi_to_ensembl, ensembl_to_ncbi):
    enhanced_entries = []
    for formbioid, proteinid, source_type, taxid, geneid in matched_entries:
        if geneid == 'None':
            enhanced_entries.append((formbioid, proteinid, source_type, taxid, None, None))
        elif source_type == 'NCBI':
            ensembl_geneid = ncbi_to_ensembl.get(geneid, None)
            enhanced_entries.append((formbioid, proteinid, source_type, taxid, geneid, ensembl_geneid))
        elif source_type == 'Ensembl':
            ncbi_geneid = ensembl_to_ncbi.get(geneid, None)
            enhanced_entries.append((formbioid, proteinid, source_type, taxid, ncbi_geneid, geneid))
        else:
            # Handle unexpected source type if necessary
            enhanced_entries.append((formbioid, proteinid, source_type, taxid, geneid, None))
    return enhanced_entries


# In[10]:


ncbi_to_ensembl, ensembl_to_ncbi = parse_geneid_mapping(mapping_file_path)


# In[11]:


enhanced_matched_entries = enhance_matched_entries(matched_entries, ncbi_to_ensembl, ensembl_to_ncbi)


# In[12]:
print("parsed ncbi to ensembl")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
def update_gene_list_with_exon_info(gff3_file_path, gene_list):
    # Step 1: Parse the GFF3 file to build a dictionary of gene IDs and their exon counts
    exon_counts = {}
    with open(gff3_file_path, 'r') as gff3_file:
        for line in gff3_file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip headers and empty lines
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                attributes = parts[8]
                # Extract gene_id from attributes
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('gene_id='):
                        gene_id = attr.split('=')[1]
                        break
                if gene_id:
                    exon_counts[gene_id] = exon_counts.get(gene_id, 0) + 1

    # Step 2: Check each gene in your list for single/multi-exon status
    updated_gene_list = []
    for gene_tuple in gene_list:
        formbioid = gene_tuple[0]
        # Assume multi-exon by default, change to False if exon count is 1 or gene is not found (implying no exons listed)
        is_multi_exon = exon_counts.get(formbioid, 0) > 1
        updated_gene_list.append(gene_tuple + (is_multi_exon,))

    return updated_gene_list


# In[13]:


updated_gene_list = update_gene_list_with_exon_info(annotation_file, enhanced_matched_entries)


# In[14]:


def parse_ensembl_orthologs(file_path):
    """
    Parses the Ensembl ortholog file to create two mappings:
    1. gene_id_to_family_id: Maps each Ensembl GeneID to its Ortholog Family ID
    2. family_id_to_gene_ids: Maps each Ortholog Family ID to a list of associated GeneIDs
    """
    gene_id_to_family_id = {}
    family_id_to_gene_ids = {}
    
    with open(file_path, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            parts = line.strip().split(',')
            family_id, gene_id = parts[0], parts[2]
            
            # Update gene_id_to_family_id mapping
            gene_id_to_family_id[gene_id] = family_id
            
            # Update family_id_to_gene_ids mapping
            if family_id not in family_id_to_gene_ids:
                family_id_to_gene_ids[family_id] = []
            family_id_to_gene_ids[family_id].append(gene_id)
    
    return gene_id_to_family_id, family_id_to_gene_ids


# In[15]:
gene_id_to_family_id, family_id_to_gene_ids = parse_ensembl_orthologs(ensembl_ortholog_file_path)


# In[16]:


def update_gene_list_with_orthologs(gene_list, gene_id_to_family_id, family_id_to_gene_ids):
    """
    Updates the gene list with a new column for Ensembl orthologs.
    """
    updated_gene_list = []
    for item in gene_list:
        ensembl_gene_id = item[5]  # Assuming this is the position of the Ensembl GeneID in the tuple
        ortholog_family_id = gene_id_to_family_id.get(ensembl_gene_id)
        
        if ortholog_family_id and ensembl_gene_id in family_id_to_gene_ids.get(ortholog_family_id, []):
            # Get all gene IDs in the family excluding the current gene ID
            orthologs = [gid for gid in family_id_to_gene_ids[ortholog_family_id]]
            updated_gene_list.append(item + (' '.join(orthologs),))
        else:
            updated_gene_list.append(item + (None,))
    
    return updated_gene_list


# In[17]:


updated_gene_list = update_gene_list_with_orthologs(updated_gene_list, gene_id_to_family_id, family_id_to_gene_ids)


# In[18]:


def parse_ncbi_orthologs(file_path):
    """
    Parses the NCBI ortholog file to create mappings similar to Ensembl:
    1. gene_id_to_family_id: Maps each NCBI GeneID and Other_GeneID to the primary GeneID (as the family ID).
    2. family_id_to_gene_ids: Maps each primary GeneID (as the family ID) to a list of associated Other_GeneIDs.
    """
    gene_id_to_family_id = {}
    family_id_to_gene_ids = {}
    
    with open(file_path, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                primary_gene_id, other_gene_id = parts[1], parts[4]
                
                # Map primary and other gene IDs to the primary gene ID as family_id
                gene_id_to_family_id[primary_gene_id] = primary_gene_id
                gene_id_to_family_id[other_gene_id] = primary_gene_id
                
                # Update family_id_to_gene_ids mapping
                if primary_gene_id not in family_id_to_gene_ids:
                    family_id_to_gene_ids[primary_gene_id] = set()
                family_id_to_gene_ids[primary_gene_id].add(other_gene_id)
                family_id_to_gene_ids[primary_gene_id].add(primary_gene_id)  # Include primary gene ID in its own group
    
    # Convert sets to lists for consistency with your requirements
    for family_id in family_id_to_gene_ids:
        family_id_to_gene_ids[family_id] = list(family_id_to_gene_ids[family_id])
    
    return gene_id_to_family_id, family_id_to_gene_ids


# In[19]:

ncbi_gene_id_to_family_id, ncbi_family_id_to_gene_ids = parse_ncbi_orthologs(ncbi_ortholog_file_path)


# In[20]:


def update_gene_list_with_ncbi_orthologs(gene_list, gene_id_to_family_id, family_id_to_gene_ids):
    updated_list = []
    for gene_tuple in gene_list:
        ncbi_gene_id = gene_tuple[4]  # Assuming NCBI GeneID is at index 4
        # Check if the NCBI GeneID is in the gene_id_to_family_id mapping
        if ncbi_gene_id in gene_id_to_family_id:
            # Find the family ID for this NCBI GeneID
            family_id = gene_id_to_family_id[ncbi_gene_id]
            # Fetch all gene IDs associated with this family ID
            gene_ids = family_id_to_gene_ids.get(family_id, [])
            # Convert the list of gene IDs into a string separated by spaces
            orthologs_str = ' '.join(gene_ids)
            updated_tuple = gene_tuple + (orthologs_str,)
        else:
            # Append None if no family ID is found for the NCBI GeneID
            updated_tuple = gene_tuple + (None,)
        
        updated_list.append(updated_tuple)
    
    return updated_list


# In[21]:


updated_gene_list_with_ncbi = update_gene_list_with_ncbi_orthologs(updated_gene_list, ncbi_gene_id_to_family_id, ncbi_family_id_to_gene_ids)


# In[22]:


def parse_ensembl_exon_counts(file_path):
    geneid_to_exon_count = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            gene_id, exon_count = parts[0], parts[2]
            geneid_to_exon_count[gene_id] = exon_count
            
    return geneid_to_exon_count


# In[23]:


ensembl_geneid_to_exon_count = parse_ensembl_exon_counts(gene2ensembl_mapping_file)


# In[24]:


def update_gene_list_with_exon_counts(updated_gene_list_with_ncbi, geneid_to_exon_count):
    updated_list = []
    for gene_tuple in updated_gene_list_with_ncbi:
        ensembl_gene_id = gene_tuple[5]  # Ensembl GeneID
        single_exon = gene_tuple[6]  # Single exon
        ensembl_orthologs = gene_tuple[7]  # All Ensembl orthologs separated by a space
        
        # Check if Ensembl GeneID, Single Exon, and Ensembl Orthologs are not 'None'
        if ensembl_gene_id != None:
            if ensembl_orthologs != None:
                # Create a set of unique gene IDs including the Ensembl GeneID and its orthologs
                gene_ids_set = set([ensembl_gene_id] + ensembl_orthologs.split(' '))
            else:
                gene_ids_set = set([ensembl_gene_id])
            
            # Counters for single and multi exon genes
            single_exon_count = 0
            multi_exon_count = 0
            
            # Iterate through the set and count based on exon numbers
            for gene_id in gene_ids_set:
                exon_count = geneid_to_exon_count.get(gene_id)
                if exon_count == '1':
                    single_exon_count += 1
                elif exon_count and int(exon_count) > 1:
                    multi_exon_count += 1
            
            # Append counts to the tuple
            updated_tuple = gene_tuple + (single_exon_count, multi_exon_count)
        else:
            # If any required data is 'None', append zeros for counts
            updated_tuple = gene_tuple + (0, 0)
        
        updated_list.append(updated_tuple)
    
    return updated_list


# In[25]:


updated_gene_list_with_ncbi_ensemblexons = update_gene_list_with_exon_counts(updated_gene_list_with_ncbi, ensembl_geneid_to_exon_count)


# In[26]:


# Function to parse the NCBI exon count file
def parse_ncbi_exon_counts(file_path):
    geneid_to_exon_count = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split(',')
            gene_id, exon_count = parts[0].strip(), parts[1].strip()
            geneid_to_exon_count[gene_id] = exon_count
    return geneid_to_exon_count


# In[27]:
ncbi_geneid_to_exon_count = parse_ncbi_exon_counts(ncbi_exon_mapping_file)

print("parsed exon counts")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
# In[28]:


# Function to update the gene list with NCBI exon counts
def update_with_ncbi_exon_counts(gene_list, exon_counts):
    updated_list = []
    for entry in gene_list:
        formbioid, proteinid, source, taxid, ncbi_geneid, ensembl_geneid, single_exon, ensembl_orthologs, ncbi_orthologs, ens_single_exon_count, ens_multi_exon_count = entry
        
        # Initialize counts
        ncbi_single_exon_count = 0
        ncbi_multi_exon_count = 0
        
        # Check if NCBI GeneID and orthologs are not None
        if ncbi_geneid != None:
            if ncbi_orthologs != None:
                gene_ids = [ncbi_geneid] + ncbi_orthologs.split()
            else:
                gene_ids = [ncbi_geneid]
            for gene_id in gene_ids:
                exon_count = exon_counts.get(gene_id)
                if exon_count:
                    if int(exon_count) == 1:
                        ncbi_single_exon_count += 1
                    else:
                        ncbi_multi_exon_count += 1
        
        # Append counts to the tuple
        updated_entry = entry + (ncbi_single_exon_count, ncbi_multi_exon_count)
        updated_list.append(updated_entry)
    
    return updated_list


# In[29]:


# Update the gene list with NCBI exon counts
updated_gene_list_with_ncbi_ensemblexons_ncbiexons = update_with_ncbi_exon_counts(updated_gene_list_with_ncbi_ensemblexons, ncbi_geneid_to_exon_count)


# # DONE PARSING ORTHOLOGS

# In[30]:
print("done with database")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
def extract_gene_ids(annotation_file):
    gene_ids = set()
    with open(annotation_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or len(row) < 9:
                continue
            feature_type = row[2]  # This is where the feature type (e.g., gene, mRNA, exon) is specified.
            if feature_type != "gene":  # Skip rows that are not genes.
                continue
            attributes = dict(field.split('=') for field in row[8].split(';') if '=' in field)
            gene_id = attributes.get('ID')
            if gene_id:
                gene_ids.add(gene_id)
    return gene_ids

def parse_gff3_for_gene_exons(gff3_file, gene_id):
    """
    Parses a GFF3 file for a specific gene ID to find its coordinates, strand, sequenceid,
    and its exon coordinates.
    """
    print(f"gene_id: {gene_id}")
    exons = []
    with open(gff3_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or len(row) < 9:
                continue
            attributes = {key: value for key, value in [field.split('=') for field in row[8].split(';') if '=' in field]}
            if 'gene_id' in attributes and attributes['gene_id'] == gene_id and row[2] == 'exon':
                sequenceid = row[0]  # Extract the sequenceid
                exons.append((sequenceid, int(row[3]), int(row[4]), row[6]))  # sequenceid, Start, End, Strand
    if exons:
        # Assuming all exons belong to the same gene and thus have the same sequenceid and strand
        return exons[0][0], exons, exons[0][3]  # sequenceid, exons, strand
    else:
        return None, [], None


# In[31]:


def parse_busco_tsv(busco_file):
    """
    Parses a BUSCO TSV file and returns a dictionary with gene IDs as keys and
    their BUSCO status (Complete or Fragmented) as values.
    """
    busco_status = {}
    with open(busco_file, 'r') as file:
        for line in file:
            if line.startswith('#'):  # Skip header or comment lines
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:  # Ensure there are enough columns
                busco_id, status, gene_id = parts[:3]
                busco_status[gene_id] = status
    return busco_status

def parse_splign_gff3(alignment_file):
    """
    Parses a Splign GFF3 file to extract alignment blocks with detailed attributes,
    including the sequence ID.
    """
    alignments = []
    with open(alignment_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or len(row) < 9:
                continue
            sequenceid = row[0]  # Extracting sequenceid
            attributes = {attr.split('=')[0]: attr.split('=')[1] for attr in row[8].split(';') if '=' in attr}
            start, end, strand = int(row[3]), int(row[4]), row[6]
            alignments.append({
                'sequenceid': sequenceid,
                'start': start,
                'end': end,
                'strand': strand
            })
    return alignments


# In[32]:


def calculate_coverage_splign(gene_sequenceid, gene_start, gene_end, alignments):
    """
    Calculates the coverage of the gene by alignment blocks and returns the coverage percentage.
    Additionally, prints the alignment lines that match the gene region, ensuring sequenceid matches.
    """
    gene_length = gene_end - gene_start + 1
    covered_segments = []
    matching_alignments = []  # To store matching alignment details

    for alignment in alignments:
        # Check for sequenceid match before assessing overlap
        if alignment['sequenceid'] == gene_sequenceid:
            overlap_start = max(alignment['start'], gene_start)
            overlap_end = min(alignment['end'], gene_end)
            if overlap_start <= overlap_end:
                covered_segments.append((overlap_start, overlap_end))
                matching_alignments.append(alignment)  # Add the matching alignment to the list

    # After processing all alignments, print the matching ones
    if matching_alignments:
        print("Matching alignment lines:")
        for match in matching_alignments:
            print(f"Sequence ID: {match['sequenceid']}, Start: {match['start']}, End: {match['end']}, Strand: {match['strand']}")

    covered_segments.sort()
    merged_segments = []
    for start, end in covered_segments:
        if not merged_segments or merged_segments[-1][1] < start:
            merged_segments.append([start, end])
        else:
            merged_segments[-1][1] = max(merged_segments[-1][1], end)

    covered_length = sum(end - start + 1 for start, end in merged_segments)
    coverage_percentage = (covered_length / gene_length) * 100
    return coverage_percentage

def calculate_coverage_splign_exons(exons, alignments):
    """
    Calculates the coverage of the gene's exons by alignment blocks and returns the coverage percentage.
    Additionally, prints the alignment lines that match the exon regions, ensuring sequenceid matches.
    """
    total_exon_length = sum(end - start + 1 for _, start, end, _ in exons)
    covered_segments = []
    matching_alignments = []
    print(exons)
    print(f"total_exon_length: {total_exon_length}")
    
    if total_exon_length > 0:
        for gene_sequenceid, exon_start, exon_end, strand in exons:
            for alignment in alignments:
                if alignment['sequenceid'] == gene_sequenceid:
                    overlap_start = max(alignment['start'], exon_start)
                    overlap_end = min(alignment['end'], exon_end)
                    if overlap_start <= overlap_end:
                        covered_segments.append((overlap_start, overlap_end))
                        matching_alignments.append(alignment)

        if matching_alignments:
            print("Matching alignment lines:")
            for match in matching_alignments:
                print(f"Sequence ID: {match['sequenceid']}, Start: {match['start']}, End: {match['end']}, Strand: {match['strand']}")

        covered_segments.sort()
        merged_segments = []
        for start, end in covered_segments:
            if not merged_segments or merged_segments[-1][1] < start:
                merged_segments.append([start, end])
            else:
                merged_segments[-1][1] = max(merged_segments[-1][1], end)

        covered_length = sum(end - start + 1 for start, end in merged_segments)
        coverage_percentage = (covered_length / total_exon_length) * 100
        return coverage_percentage
    else:
        return 0


# In[33]:


def calculate_coverage_pasa_exons(exons, alignments):
    """
    Calculates the coverage of the gene's exons by alignment blocks and returns the coverage percentage.
    Additionally, prints the alignment lines that match the exon regions, ensuring sequenceid matches.
    """
    total_exon_length = sum(end - start + 1 for _, start, end, _ in exons)
    covered_segments = []
    matching_alignments = []
    print(exons)
    print(f"total_exon_length: {total_exon_length}")
    
    if total_exon_length > 0:
        for gene_sequenceid, exon_start, exon_end, strand in exons:
            for alignment in alignments:
                if alignment['sequenceid'] == gene_sequenceid:
                    overlap_start = max(alignment['start'], exon_start)
                    overlap_end = min(alignment['end'], exon_end)
                    if overlap_start <= overlap_end:
                        covered_segments.append((overlap_start, overlap_end))
                        matching_alignments.append(alignment)

        if matching_alignments:
            print("Matching alignment lines:")
            for match in matching_alignments:
                print(f"Sequence ID: {match['sequenceid']}, Start: {match['start']}, End: {match['end']}, Strand: {match['strand']}")

        covered_segments.sort()
        merged_segments = []
        for start, end in covered_segments:
            if not merged_segments or merged_segments[-1][1] < start:
                merged_segments.append([start, end])
            else:
                merged_segments[-1][1] = max(merged_segments[-1][1], end)

        covered_length = sum(end - start + 1 for start, end in merged_segments)
        coverage_percentage = (covered_length / total_exon_length) * 100
        return coverage_percentage
    else:
        return 0


# In[34]:


gene_ids = extract_gene_ids(annotation_file)


# In[35]:


# In[36]:


splign_alignments = parse_splign_gff3(splign_file) if splign_file else []
reference_splign_alignments = parse_splign_gff3(reference_splign_file) if reference_splign_file else []
pasa_alignments = parse_splign_gff3(pasa_file) if pasa_file else []
reference_pasa_alignments = parse_splign_gff3(reference_pasa_file) if reference_pasa_file else []
busco_status = parse_busco_tsv(busco_file) if busco_file else {}

print("parsed alignment files")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
# In[37]:


def update_gene_list_with_annotation_exons(annotation_file_path, gene_list):
    # Parse the GFF3 file to build a dictionary of formbioid and their exon counts
    formbioid_to_exon_counts = {}
    with open(annotation_file_path, 'r') as gff3_file:
        for line in gff3_file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip headers and empty lines
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                attributes = parts[8]
                formbioid = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        formbioid = attr.split('=')[1]
                        break
                if formbioid:
                    formbioid_to_exon_counts[formbioid] = formbioid_to_exon_counts.get(formbioid, 0) + 1

    # Extract formbioids from the gene list
    existing_formbioids = set([gene_tuple[0] for gene_tuple in gene_list])
    added_list = []
    # Update the gene list by adding missing formbioids with their multi-exon status
    for formbioid, exon_count in formbioid_to_exon_counts.items():
        if formbioid not in existing_formbioids:
            print(f"{formbioid} {exon_count}")
            # Determine multi-exon status
            is_multi_exon = exon_count > 1
            # Create a new tuple with null for missing fields except formbioid and multi-exon status
            new_tuple = (formbioid, None, None, None, None, None, is_multi_exon, None, None, 0, 0, 0, 0)
            added_list.append((formbioid, exon_count))
            gene_list.append(new_tuple)

    return gene_list, added_list


# In[38]:


def parse_annotation_for_exon_counts(annotation_file_path):
    """
    Parses the GFF3 annotation file to extract formbioids and determine their exon counts.
    
    Parameters:
    - annotation_file_path: Path to the GFF3 annotation file.
    
    Returns:
    A dictionary mapping formbioids to their exon counts.
    """
    formbioid_to_exon_counts = {}
    with open(annotation_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Skip headers and empty lines
            parts = line.strip().split('\t')
            if parts[2] == 'exon':
                attributes = parts[8]
                formbioid = None
                for attr in attributes.split(';'):
                    if attr.startswith('Parent='):
                        formbioid = attr.split('=')[1].split('.')[0]  # Splitting to remove transcript identifier
                        # If parent_id contains '.t', strip it to get the gene formbioid
                        if ".t" in formbioid:
                            formbioid = formbioid.split('.t')[0]
                        break
                if formbioid:
                    formbioid_to_exon_counts[formbioid] = formbioid_to_exon_counts.get(formbioid, 0) + 1
    
    return formbioid_to_exon_counts

def update_gene_list_with_missing_formbioids(gene_list, formbioid_to_exon_counts):
    """
    Updates the gene list with any missing formbioids and their multi-exon status based on exon counts.
    
    Parameters:
    - gene_list: The existing list of gene tuples.
    - formbioid_to_exon_counts: A dictionary mapping formbioids to their exon counts.
    
    Returns:
    An updated list of gene tuples including missing formbioids with their multi-exon status.
    """
    existing_formbioids = {gene_tuple[0] for gene_tuple in gene_list}
    updated_list = gene_list[:]
    added_list = []
    for formbioid, exon_count in formbioid_to_exon_counts.items():
        if formbioid not in existing_formbioids:
            # Determine multi-exon status based on exon count
            is_multi_exon = exon_count > 1
            new_tuple = (formbioid, None, None, None, None, None, is_multi_exon, None, None, 0, 0, 0, 0)
            updated_list.append(new_tuple)
            added_list.append((formbioid, exon_count))
    
    return updated_list, added_list


# In[39]:


formbioid_to_multi_exon_counts = parse_annotation_for_exon_counts(annotation_file)


# In[40]:


updated_gene_list_with_missing_formbioids, added_list = update_gene_list_with_missing_formbioids(updated_gene_list_with_ncbi_ensemblexons_ncbiexons, formbioid_to_multi_exon_counts)

print("updated with missing ids")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
# # Add buscos

# In[41]:


genes_in_busco = {gene_id for gene_id, status in busco_status.items() if status in ["Complete", "Fragmented"]}
# Use a set comprehension to remove '.t1' from the end of each gene ID where present
genes_in_busco = {gene_id[:-3] if gene_id.endswith('.t1') else gene_id for gene_id in genes_in_busco}


# In[42]:


updated_gene_list_with_busco_status = [
    entry + (entry[0] in genes_in_busco,) for entry in updated_gene_list_with_missing_formbioids
]


# # Checking pasa alignments

# In[43]:


def build_gene_to_exons_map(gff3_file):
    gene_to_exons = {}
    with open(gff3_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or len(row) < 9:
                continue
            attributes = dict(field.split('=') for field in row[8].split(';') if '=' in field)
            if 'gene_id' in attributes and row[2] == 'exon':
                gene_id = attributes['gene_id']
                exon_info = (row[0], int(row[3]), int(row[4]), row[6])  # sequenceid, Start, End, Strand
                if gene_id not in gene_to_exons:
                    gene_to_exons[gene_id] = []
                gene_to_exons[gene_id].append(exon_info)
    return gene_to_exons

# Use the map for fast lookup
gene_to_exons_map = build_gene_to_exons_map(annotation_file)


# In[45]:


def calculate_coverage_splign_exons_from_map(gene_id, gene_to_exons_map, alignments):
    """
    Calculates the coverage of a gene's exons by alignment blocks using a pre-processed map
    from gene IDs to exon information, returning the coverage percentage.
    """
    if gene_id not in gene_to_exons_map:
        return 0  # or handle missing gene_id appropriately

    exons = gene_to_exons_map[gene_id]
    total_exon_length = sum(exon[2] - exon[1] + 1 for exon in exons)
    covered_segments = []
    
    for exon in exons:
        gene_sequenceid, exon_start, exon_end, strand = exon
        for alignment in alignments:
            if alignment['sequenceid'] == gene_sequenceid:
                overlap_start = max(alignment['start'], exon_start)
                overlap_end = min(alignment['end'], exon_end)
                if overlap_start <= overlap_end:
                    covered_segments.append((overlap_start, overlap_end))

    # Merge overlapping or contiguous segments
    covered_segments.sort(key=lambda x: x[0])
    merged_segments = []
    for start, end in covered_segments:
        if not merged_segments or merged_segments[-1][1] < start - 1:
            merged_segments.append([start, end])
        else:
            merged_segments[-1][1] = max(merged_segments[-1][1], end)

    # Calculate covered length from merged segments
    covered_length = sum(end - start + 1 for start, end in merged_segments)
    coverage_percentage = (covered_length / total_exon_length) * 100 if total_exon_length > 0 else 0

    return coverage_percentage


# In[46]:


def update_gene_list_with_alignment_coverage(updated_gene_list, gene_to_exons_map, alignments):
    updated_list = []
    for gene_tuple in updated_gene_list:
        gene_id = gene_tuple[0]
        coverage_percentage = calculate_coverage_splign_exons_from_map(gene_id, gene_to_exons_map, alignments)
        updated_list.append(gene_tuple + (coverage_percentage,))
    return updated_list


# In[47]:


def index_alignments_by_sequenceid(alignments):
    indexed_alignments = {}
    for alignment in alignments:
        sequenceid = alignment['sequenceid']
        if sequenceid not in indexed_alignments:
            indexed_alignments[sequenceid] = []
        indexed_alignments[sequenceid].append(alignment)
    return indexed_alignments


# In[48]:


# Index alignments by sequenceid for faster access
indexed_alignments = index_alignments_by_sequenceid(splign_alignments)


# In[49]:


def calculate_coverage_splign_exons_from_map(gene_id, gene_to_exons_map, indexed_alignments):
    if gene_id not in gene_to_exons_map:
        return 0  # Gene ID not found in the map

    exons = gene_to_exons_map[gene_id]
    total_exon_length = sum(exon_end - exon_start + 1 for _, exon_start, exon_end, _ in exons)
    covered_segments = []

    for exon in exons:
        gene_sequenceid, exon_start, exon_end, _ = exon
        for alignment in indexed_alignments.get(gene_sequenceid, []):
            overlap_start = max(alignment['start'], exon_start)
            overlap_end = min(alignment['end'], exon_end)
            if overlap_start <= overlap_end:
                covered_segments.append((overlap_start, overlap_end))

    # Merge overlapping segments
    covered_segments.sort(key=lambda x: x[0])
    merged_segments = []
    for start, end in covered_segments:
        if not merged_segments or merged_segments[-1][1] < start - 1:
            merged_segments.append([start, end])
        else:
            # Extend the end of the last segment if the current segment overlaps or is consecutive
            merged_segments[-1][1] = max(merged_segments[-1][1], end)

    # Calculate total covered length from merged segments
    covered_length = sum(end - start + 1 for start, end in merged_segments)
    coverage_percentage = (covered_length / total_exon_length) * 100 if total_exon_length > 0 else 0
    return min(coverage_percentage, 100)  # Ensure coverage does not exceed 100%

def update_gene_list_with_alignment_coverage(gene_list, gene_to_exons_map, indexed_alignments):
    updated_list = []
    for gene_tuple in gene_list:
        gene_id = gene_tuple[0]
        coverage_percentage = calculate_coverage_splign_exons_from_map(gene_id, gene_to_exons_map, indexed_alignments)
        updated_list.append(gene_tuple + (coverage_percentage,))
    return updated_list

# Use the indexed alignments for updating the gene list
updated_gene_list_with_coverage = update_gene_list_with_alignment_coverage(
    updated_gene_list_with_busco_status, gene_to_exons_map, indexed_alignments)


# # Reference splign

# In[50]:


indexed_alignments = index_alignments_by_sequenceid(reference_splign_alignments)


# In[51]:


# Use the indexed alignments for updating the gene list
updated_gene_list_with_coverage_allsplign = update_gene_list_with_alignment_coverage(
    updated_gene_list_with_coverage, gene_to_exons_map, indexed_alignments)


# # Now for pasa and reference pasa

# In[52]:


indexed_alignments = index_alignments_by_sequenceid(pasa_alignments)


# In[53]:


updated_gene_list_with_coverage_allsplign_pasa = update_gene_list_with_alignment_coverage(
    updated_gene_list_with_coverage_allsplign, gene_to_exons_map, indexed_alignments)


# In[54]:


indexed_alignments = index_alignments_by_sequenceid(reference_pasa_alignments)


# In[55]:


updated_gene_list_with_coverage_allsplign_allpasa = update_gene_list_with_alignment_coverage(
    updated_gene_list_with_coverage_allsplign_pasa, gene_to_exons_map, indexed_alignments)


# # Now for miniprot

# In[56]:


def parse_protein_alignment_exons(protein_alignment_file):
    protein_exons = {}
    with open(protein_alignment_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            feature_type = parts[2]
            if feature_type == 'exon':  # Assuming you're interested in both exons and CDS
                sequenceid = parts[0]
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                # Parse attributes
                # attributes = {attr.strip().split('=')[0]: attr.strip().split('=')[1].replace('"', '') for attr in parts[8].split(';') if '=' in attr}
                # gene_id = attributes.get('gene_id')

                # Store exon information in a dictionary with lists of dictionaries for each sequenceid
                if sequenceid not in protein_exons:
                    protein_exons[sequenceid] = []
                protein_exons[sequenceid].append({
                    'sequenceid': sequenceid,
                    'start': start,
                    'end': end,
                    'strand': strand
                })
    return protein_exons


# In[57]:
if protein_alignment_file != None:
    protein_exons = parse_protein_alignment_exons(protein_alignment_file)
else:
    protein_exons = {}


# In[58]:


finalTuples = update_gene_list_with_alignment_coverage(
    updated_gene_list_with_coverage_allsplign_allpasa, gene_to_exons_map, protein_exons)


# In[59]:


def create_csv_from_tuples(tuples_list, csv_file_path):
    """
    Creates a CSV file from a list of tuples with headers.
    
    Parameters:
    - tuples_list: List of tuples to write into the CSV.
    - csv_file_path: Path to save the CSV file.
    """
    headers = [
        'Formbioid', 'Protein Accession', 'Eggnog Ortholog Source', 'TaxID', 'NCBI GeneID', 'Ensembl GeneID',
        'Multi Exon T/F', 'Ensembl Orthologs', 'NCBI Orthologs', 'Ensembl Single Exon Count',
        'Ensembl Multi Exon Count', 'NCBI Single Exon Count', 'NCBI Multi Exon Count',
        'Busco T/F', 'Splign Alignment Percentage', 'Reference Splign Alignment Percentage', 'Pass Alignment Percentage',
        'Reference Pass Alignment Percentage', 'Miniprot Alignment Percentage'
    ]

    # Sort the list by formbioid
    sorted_list = sorted(tuples_list, key=lambda x: x[0])
    
    # Write the sorted list to a CSV file
    with open(csv_file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(sorted_list)

    return csv_file_path


# In[60]:


def add_reasons_to_full_list(full_list, multi_exon_conditions, single_exon_conditions):
    """
    Adds a reason the full_list based on provided conditions.

    Parameters:
    - full_list: List of tuples with the structure [formbioid, protein accession, source type, taxid, ncbi geneid, ensembl geneid, Multi exon true or false, Ensembl orthologs, ncbi orthologs, ensembl single exon count, ensembl multi exon count, ncbi single exon count, ncbi multi exon count, isBusco True or False, Splign alignment, reference Splign alignment, pass alignment, reference pass alignment, miniprot]
    - conditions: Dictionary with keys as condition types and values as the condition. Example: {"is_single_exon": True, "has_multi_exon_orthologs": True, "splign_coverage_min": 90, "is_busco": True}

    Returns:
    An updated list of tuples with reasons added at the end of each tuple.
    """
    updated_list = []
    for item in full_list:
        reasons = []
        # Extract relevant fields for easy comparison
        ortholog_found, is_multi_exon, ensembl_orthologs, ncbi_orthologs, ensembl_single_exon_count, ensembl_multi_exon_count, ncbi_single_exon_count, ncbi_multi_exon_count, is_busco, splign_alignment, reference_splign_alignment, pasa_alignment, reference_pasa_alignment, miniprot_alignment = item[1], item[6], item[7], item[8], item[9], item[10], item[11], item[12], item[13], item[14], item[15], item[16], item[17], item[18]
        total_single_exon = ensembl_single_exon_count + ncbi_single_exon_count
        total_multi_exon = ensembl_multi_exon_count + ncbi_multi_exon_count
        total_orthologs = total_single_exon + total_multi_exon
        if is_multi_exon:
            # Gene is multi exon so should use multi exon conditions
            if total_orthologs != 0:
                percent_same_exon_count = (total_multi_exon / total_orthologs) * 100
            else:
                percent_same_exon_count = 0
            
            conditions = multi_exon_conditions
            
        else:
            # Gene is single exon so should use single exon conditions
            if total_orthologs != 0:
                percent_same_exon_count = (total_single_exon / total_orthologs) * 100
            else:
                percent_same_exon_count = 0
                
            conditions = single_exon_conditions
            # print(f"{ortholog_found}: single")
           
        if "num_orthologs" in conditions:
            num_orthologs_found = 0
            if ensembl_orthologs == None and ncbi_orthologs == None:
                if ortholog_found == None:
                    num_orthologs_found = 0
                else:
                    num_orthologs_found = 1
            else:
                if ensembl_orthologs != None:
                    num_orthologs_found += len(ensembl_orthologs.split())
                if ncbi_orthologs != None:
                    num_orthologs_found += len(ncbi_orthologs.split())

            if conditions['num_orthologs'] <= num_orthologs_found:
                if "percent_same_exon_match" in conditions:
                    if conditions['percent_same_exon_match'] <= percent_same_exon_count:
                        if is_multi_exon:
                            # print(f"multi: {conditions['percent_same_exon_match']} <= {percent_same_exon_count} and {ensembl_single_exon_count}, {ensembl_multi_exon_count}, {ncbi_single_exon_count}, {ncbi_multi_exon_count}")
                            reasons.append(f"{num_orthologs_found} orthologs found")
                            reasons.append(f"{percent_same_exon_count}% of orthologs are also multi exon")
                        else:
                            # print(f"{conditions['percent_same_exon_match']} <= {percent_same_exon_count} and {ensembl_single_exon_count}, {ensembl_multi_exon_count}, {ncbi_single_exon_count}, {ncbi_multi_exon_count}")
                            reasons.append(f"{num_orthologs_found} orthologs found")
                            reasons.append(f"{percent_same_exon_count}% of orthologs are also single exon")
                    # else:
                    #     print(f"NO: {conditions['percent_same_exon_match']} <= {percent_same_exon_count} and {ensembl_single_exon_count}, {ensembl_multi_exon_count}, {ncbi_single_exon_count}, {ncbi_multi_exon_count}")
                else:
                    reasons.append(f"{num_orthologs_found} orthologs found")
        else:
            if "percent_same_exon_match" in conditions and conditions['percent_same_exon_match'] <= percent_same_exon_count:
                if is_multi_exon:
                    reasons.append(f"{percent_same_exon_count}% of orthologs are also multi exon")
                else:
                    # print(f"{conditions['percent_same_exon_match']} <= {percent_same_exon_count} and {ensembl_single_exon_count}, {ensembl_multi_exon_count}, {ncbi_single_exon_count}, {ncbi_multi_exon_count}")
                    reasons.append(f"{percent_same_exon_count}% of orthologs are also single exon")
        
        if "is_busco" in conditions and bool(conditions['is_busco']) == True and bool(is_busco) == True:
            reasons.append(f"Busco match")
            
        if "splign_coverage_min" in conditions and conditions['splign_coverage_min'] <= splign_alignment:
            reasons.append(f"{splign_alignment}% of exons are covered by splign transcript alignments")
        
        if "splign_coverage_min" in conditions and conditions['reference_splign_coverage_min'] <= reference_splign_alignment:
            reasons.append(f"{reference_splign_alignment}% of exons are covered by reference splign transcript alignments")
            
        if "pasa_coverage_min" in conditions and conditions['pasa_coverage_min'] <= pasa_alignment:
            reasons.append(f"{pasa_alignment}% of exons are covered by pasa transcript alignments")
        
        if "reference_pasa_coverage_min" in conditions and conditions['reference_pasa_coverage_min'] <= reference_pasa_alignment:
            reasons.append(f"{reference_pasa_alignment}% of exons are covered by reference pasa transcript alignments")
            
        if "miniprot_coverage_min" in conditions and conditions['miniprot_coverage_min'] <= miniprot_alignment:
            reasons.append(f"{miniprot_alignment}% of exons are covered by miniprot protein alignments")
            # if is_multi_exon == False:
            #     print(f"{conditions['miniprot_coverage_min']} <= {miniprot_alignment}")
        
        # if "reference_pasa_coverage_min" in conditions and conditions['reference_pasa_coverage_min'] <= reference_pasa_alignment:
        #     reason = f"{reason} {reference_pasa_alignment}% of exons are covered by reference pasa transcript alignments,"

        # Add more conditions as needed
        
        # Remove any white space at the beginning
        # reason = reason.lstrip()

        # # Remove trailing commas
        # reason = reason.rstrip(',')

        # Add the reason string to the end of the tuple
        updated_item = item + ("| ".join(reasons),)
        updated_list.append(updated_item)
        
    return updated_list


# In[61]:
def adjust_conditions_with_logic(single_exon_conditions, adjustment_levels):
    """
    Adjusts single_exon_conditions with specific logic for is_busco and num_orthologs.

    Parameters:
    - single_exon_conditions: The conditions to be adjusted.
    - adjustment_levels: A dictionary specifying how much to adjust each condition.
    """

    # Flag to indicate if splign_coverage_min reaches 0
    adjust_num_orthologs = False
    add_miniprot = False
    # Adjust conditions based on adjustment levels except for is_busco
    for condition, adjustment in adjustment_levels.items():
        if condition == 'is_busco':
            continue  # Skip adjusting is_busco, it should always remain True
        if condition in single_exon_conditions and isinstance(single_exon_conditions[condition], int):
            new_value = max(0, single_exon_conditions[condition] + adjustment)
   
            if condition == 'splign_coverage_min' and new_value <= 0:
                adjust_num_orthologs = True  # Enable num_orthologs adjustment if splign_coverage_min reaches 0
            if condition == 'splign_coverage_min' and new_value <= 25:
                add_miniprot = True
            if (condition == 'splign_coverage_min' or condition == 'pasa_coverage_min' or condition == 'reference_pasa_coverage_min' or condition == 'reference_splign_coverage_min') and new_value <= 0:
                if single_exon_conditions['miniprot_coverage_min'] > 0:
                    new_value = 1
            single_exon_conditions[condition] = new_value
            
    if add_miniprot:
        # Ensure 'miniprot_coverage_min' exists, add it with a default value if not
        if 'miniprot_coverage_min' not in single_exon_conditions:
            single_exon_conditions['miniprot_coverage_min'] = 100  # Default value, adjust as necessary

    # Adjust num_orthologs only if splign_coverage_min has reached 0
    if adjust_num_orthologs:
        # Only adjust num_orthologs to 0 if miniprot_coverage_min is already 0
        if single_exon_conditions['miniprot_coverage_min'] == 0:
            single_exon_conditions['num_orthologs'] = max(0, single_exon_conditions['num_orthologs'] - 1)
        else:
            # Otherwise, ensure num_orthologs does not drop below 1
            single_exon_conditions['num_orthologs'] = max(1, single_exon_conditions['num_orthologs'] - 1)

    return single_exon_conditions


def adjust_and_update_until_condition_met(finalTuples, multi_exon_conditions, single_exon_conditions, adjustment_levels, target_ratio_single_multi):
    condition_met = False
    best_difference = float('inf')  # Initialize with a large number
    best_list = []
    best_single_exon_conditions = single_exon_conditions.copy()
    best_counts = {}  # To store the counts for the best case

    while not condition_met:
        updated_list = add_reasons_to_full_list(finalTuples, multi_exon_conditions, single_exon_conditions)

        single_exon_with_reasons = [item for item in updated_list if item[-1] and not item[6]]
        multi_exon_with_reasons = [item for item in updated_list if item[-1] and item[6]]
        total_single_exon = [item for item in updated_list if not item[6]]

        num_single_exon_with_reasons = len(single_exon_with_reasons)
        num_multi_exon_with_reasons = len(multi_exon_with_reasons)
        num_total_single_exon = len(total_single_exon)

        current_ratio = num_single_exon_with_reasons / num_multi_exon_with_reasons if num_multi_exon_with_reasons > 0 else 0
        current_difference = abs(current_ratio - target_ratio_single_multi)

        # Check if current setup is closer to the target ratio than the previous best
        if current_difference < best_difference:
            best_difference = current_difference
            best_list = updated_list.copy()
            best_single_exon_conditions = single_exon_conditions.copy()
            # Update the best counts
            best_counts = {
                'best_single_exon_with_reasons': num_single_exon_with_reasons,
                'best_multi_exon_with_reasons': num_multi_exon_with_reasons,
                'best_total_single_exon': num_total_single_exon
            }

        # Adjust conditions if necessary
        print(f"target_ratio_single_multi: {target_ratio_single_multi}")
        print(f"num_multi_exon_with_reasons: {num_multi_exon_with_reasons}")
        print("checking num_single_exon_with_reasons < target_ratio_single_multi * num_multi_exon_with_reasons")
        print(f"{num_single_exon_with_reasons} < {target_ratio_single_multi * num_multi_exon_with_reasons}")
        print("num_single_exon_with_reasons < num_total_single_exon")
        print(f"{num_single_exon_with_reasons} < {num_total_single_exon}")
        if num_single_exon_with_reasons < target_ratio_single_multi * num_multi_exon_with_reasons and num_single_exon_with_reasons < num_total_single_exon:
            print("Adjusting single exon conditions...")
            single_exon_conditions = adjust_conditions_with_logic(single_exon_conditions, adjustment_levels)
            print(single_exon_conditions)
        else:
            condition_met = True
            print("conditions met")

    return best_list, best_single_exon_conditions, best_counts




adjustment_levels = {
    "percent_same_exon_match": -1,
    "splign_coverage_min": -6,
    "reference_splign_coverage_min": -3,
    "pasa_coverage_min": -6,
    "reference_pasa_coverage_min": -3,
    "miniprot_coverage_min": -6,
    # No direct adjustment for is_busco and num_orthologs here, handled within the function
}
print("finished parsing alignment coverage")
updated_list, best_single_exon_conditions, best_counts = adjust_and_update_until_condition_met(finalTuples, multi_exon_conditions, single_exon_conditions, adjustment_levels, target_ratio_single_multi)
# Example usage
# updated_list = add_reasons_to_full_list(finalTuples, multi_exon_conditions, single_exon_conditions)

print(best_single_exon_conditions)

print(best_counts)
print("finished the update loop")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
# In[62]:


# Count the number of entries with an empty reasons column
# empty_reasons_count = sum(1 for item in updated_list if not item[-1])

# empty_reasons_count


# In[63]:


# Filter for tuples where the last column is an empty string
empty_reasons_list = [item for item in updated_list if not item[-1]]


# In[64]:


# Filter the list for tuples where the last column is an empty string and column 6 is False
# filtered_list = [item for item in updated_list if not item[-1] and item[6] == False]


# In[65]:


# Filter the list for tuples where the last column is an empty string and column 6 is False
# final_single_exon_list = [item for item in updated_list if item[6] == False]


# In[66]:


# Example usage
# Call the function
create_csv_from_tuples(updated_list, csv_file_path)


# In[67]:


headers = [
    'Formbioid', 'Protein Accession', 'Eggnog Ortholog Source', 'TaxID', 'NCBI GeneID', 'Ensembl GeneID',
    'Multi Exon T/F', 'Ensembl Orthologs', 'NCBI Orthologs', 'Ensembl Single Exon Count',
    'Ensembl Multi Exon Count', 'NCBI Single Exon Count', 'NCBI Multi Exon Count',
    'Busco T/F', 'Splign Alignment Percentage', 'Reference Splign Alignment Percentage', 'Pass Alignment Percentage',
    'Reference Pass Alignment Percentage', 'Miniprot Alignment Percentage', 'Filtered Support'
]


# In[68]:


# Filter for tuples where the last column is an empty string
supported_reasons_list = [item for item in updated_list if item[-1]]


# In[69]:


def add_model_evidence_to_gff(annotation_file, supported_list, output_file):
    # Create a dictionary for quick lookup of formbioid to Filtered Support reason
    formbioid_to_reason = {formbioid: reason for formbioid, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, reason in supported_list}
    
    with open(annotation_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                # Write out comments and whitespace lines unchanged
                outfile.write(line)
                continue
            
            parts = line.strip().split('\t')
            # Ensure there are 9 columns in a valid GFF3 line
            if len(parts) < 9:
                outfile.write(line)
                continue
            
            attributes = parts[8]
            formbioid = None
            new_attributes = attributes
            
            # Extract formbioid from the ID attribute
            for attribute in attributes.split(';'):
                if attribute.startswith('ID='):
                    formbioid = attribute.split('=')[1].split('.')[0]  # Adjust based on your ID formatting
                    break
            
            # If the formbioid is in supported_list, append the model_evidence attribute
            if formbioid and formbioid in formbioid_to_reason:
                reason = formbioid_to_reason[formbioid]
                new_attributes += f";model_evidence={reason}"
            
            # Construct the new line with updated attributes
            new_line = '\t'.join(parts[:8] + [new_attributes]) + '\n'
            outfile.write(new_line)


# In[70]:
add_model_evidence_to_gff(annotation_file,supported_reasons_list,output_evidence_gff)


# In[71]:


def write_annotations(annotation_file, removed_list, output_file):
    # Convert removed_list to a set for faster lookup
    removed_formbioids = {item[0] for item in removed_list}  # Assuming the formbioid is the first element in each tuple

    with open(annotation_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)  # Copy over GFF header and comments
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Skip incomplete lines

            attributes = parts[8]
            attribute_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}

            # Check both ID and Parent for formbioid presence
            formbioid_id = attribute_dict.get('ID', '').split('.')[0]  # Adjust as necessary based on your ID format
            formbioid_parent = attribute_dict.get('Parent', '').split('.')[0]

            # Write line to output file if either ID or Parent matches a formbioid in the removed list
            if formbioid_id in removed_formbioids or formbioid_parent in removed_formbioids:
                outfile.write(line)

# Call the function with paths to your files


# In[72]:
def filter_genes_by_length(annotation_file, supported_list, empty_list, length_threshold):
    """
    Filters genes by length, moving entries from supported_list to empty_list if they exceed the threshold.

    Parameters:
    - annotation_file: Path to the GFF3 annotation file.
    - supported_list: List of tuples representing supported genes.
    - empty_list: List of tuples representing genes without support.
    - length_threshold: Length threshold for filtering genes.

    Returns:
    - updated_supported_list: List of tuples for genes under the threshold.
    - updated_empty_list: Updated list of tuples for genes without support or exceeding length threshold.
    """
    # Parse the GFF3 file to get gene lengths
    gene_to_length = {}
    updated_empty_list = empty_list
    with open(annotation_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'gene':
                gene_id = None
                for attr in parts[8].split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                if gene_id:
                    start, end = int(parts[3]), int(parts[4])
                    gene_length = end - start + 1
                    gene_to_length[gene_id] = gene_length

    # Filter supported_list based on length
    updated_supported_list = []
    matches = 0
    for gene_tuple in supported_list:
        formbioid = gene_tuple[0]
        if formbioid in gene_to_length:
            matches += 1
            #print(gene_to_length[formbioid])
        if formbioid in gene_to_length and int(gene_to_length[formbioid]) > length_threshold:
            # Move to empty_list if gene length exceeds threshold
            updated_empty_list.append(gene_tuple)
            print("removing: gene_tuple")
        else:
            updated_supported_list.append(gene_tuple)
    # print(f"matches: {matches}")
    return updated_supported_list, updated_empty_list, gene_to_length

supported_reasons_list, empty_reasons_list, gene_to_length = filter_genes_by_length(
    annotation_file, supported_reasons_list, empty_reasons_list, length_threshold
)

write_annotations(annotation_file,empty_reasons_list,output_removed_gff)


# In[73]:
write_annotations(output_evidence_gff,supported_reasons_list,output_supported_gff)

print("Done")
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
