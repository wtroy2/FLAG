#!/bin/bash
#exonerate_p2g.sh

usage() {
  echo "-h Help documentation for tasser.sh"
  echo "-i  --input entap tsv"
  echo "-a  --input annotation gtf file"
  echo "-s  --species name"
  echo "-l  --busco lineage"
  echo "-g  --genome assembly file"
  echo "-p  --program (eggnog or entap)"
  echo "-d  --flag annotation checker database"
  echo "-b  --multi_percent_same_exon_match"
  echo "-c  --multi_splign_coverage_min"
  echo "-e  --multi_reference_splign_coverage_min"
  echo "-f  --multi_pasa_coverage_min"
  echo "-j  --multi_reference_pasa_coverage_min"
  echo "-k  --multi_miniprot_coverage_min"
  echo "-m  --multi_is_busco"
  echo "-n  --multi_num_orthologs"
  echo "-o  --single_percent_same_exon_match"
  echo "-q  --single_splign_coverage_min"
  echo "-r  --single_reference_splign_coverage_min"
  echo "-t  --single_pasa_coverage_min"
  echo "-u  --single_reference_pasa_coverage_min"
  echo "-v  --single_is_busco"
  echo "-w  --single_num_orthologs"
  echo "-y  --single_miniprot_coverage_min"
  echo "-z  --single_multi_ratio"
  echo "-x  --length_threshold"
  echo "Example: bash exonerate_p2g.sh -p GCF_000001905.1_Loxafr3.0_protein.faa -g GCF_000001905.1_Loxafr3.0_genomic.fna -t 24"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :i:a:s:l:g:p:d:b:c:e:f:j:k:m:n:o:q:r:t:u:v:w:y:z:x:h opt
do
    case $opt in
        i) entaptsv=$OPTARG;;
        a) annotation=$OPTARG;;
        s) speciesName=$OPTARG;;
        l) lineage=$OPTARG;;
        g) genome=$OPTARG;;
        p) program=$OPTARG;;
        d) flagdb=$OPTARG;;
        b) multi_percent_same_exon_match=$OPTARG;;
        c) multi_splign_coverage_min=$OPTARG;;
        e) multi_reference_splign_coverage_min=$OPTARG;;
        f) multi_pasa_coverage_min=$OPTARG;;
        j) multi_reference_pasa_coverage_min=$OPTARG;;
        k) multi_miniprot_coverage_min=$OPTARG;;
        m) multi_is_busco=$OPTARG;;
        n) multi_num_orthologs=$OPTARG;;
        o) single_percent_same_exon_match=$OPTARG;;
        q) single_splign_coverage_min=$OPTARG;;
        r) single_reference_splign_coverage_min=$OPTARG;;
        t) single_pasa_coverage_min=$OPTARG;;
        u) single_reference_pasa_coverage_min=$OPTARG;;
        v) single_is_busco=$OPTARG;;
        w) single_num_orthologs=$OPTARG;;
        y) single_miniprot_coverage_min=$OPTARG;;
        z) single_multi_ratio=$OPTARG;;
        x) length_threshold=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${entaptsv} ]] || [[ -z ${annotation} ]] || [[ -z ${genome} ]]
then
    usage
fi
if [[ -z ${speciesName} ]]
then
    speciesName=Sample_species
fi
if [[ -s ${flagdb} ]]
then
    tar -xf $flagdb
    #mkdir -p flagdb
    #tar xvfz ${flagdb} --strip-components=1 -C flagdb
    mv flagdb/* .
else
    echo "Missing Database File"
    usage
fi
if [[ -z ${threads} ]]
then
    threads=`nproc`
fi
if [[ -z ${program} ]]
then
    program="eggnog"
fi
if [[ -z ${multi_percent_same_exon_match} ]]
then
    multi_percent_same_exon_match="30"
fi
if [[ -z ${multi_splign_coverage_min} ]]
then
    multi_splign_coverage_min="75"
fi
if [[ -z ${multi_reference_splign_coverage_min} ]]
then
    multi_reference_splign_coverage_min="65"
fi
if [[ -z ${multi_pasa_coverage_min} ]]
then
    multi_pasa_coverage_min="75"
fi
if [[ -z ${multi_reference_pasa_coverage_min} ]]
then
    multi_reference_pasa_coverage_min="65"
fi
if [[ -z ${multi_miniprot_coverage_min} ]]
then
    multi_miniprot_coverage_min="85"
fi
if [[ -z ${multi_is_busco} ]]
then
    multi_is_busco="True"
fi
if [[ -z ${multi_num_orthologs} ]]
then
    multi_num_orthologs="1"
fi
if [[ -z ${single_percent_same_exon_match} ]]
then
    single_percent_same_exon_match="10"
fi
if [[ -z ${single_splign_coverage_min} ]]
then
    single_splign_coverage_min="85"
fi
if [[ -z ${single_reference_splign_coverage_min} ]]
then
    single_reference_splign_coverage_min="75"
fi
if [[ -z ${single_pasa_coverage_min} ]]
then
    single_pasa_coverage_min="85"
fi
if [[ -z ${single_reference_pasa_coverage_min} ]]
then
    single_reference_pasa_coverage_min="75"
fi
if [[ -z ${single_is_busco} ]]
then
    single_is_busco="True"
fi
if [[ -z ${single_num_orthologs} ]]
then
    single_num_orthologs="5"
fi
if [[ -z ${single_multi_ratio} ]]
then
    single_multi_ratio="0.12"
fi
if [[ -z ${length_threshold} ]]
then
    length_threshold="10000000"
fi

#checking if the fasta is a .gz
faEnd="${genome: -3}"
if [[ "${faEnd}" == ".gz" ]]; then
    cp $genome genome.fa.gz
    gunzip genome.fa.gz
else
    cp $genome genome.fa
fi

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

baseDir="`dirname \"$0\"`"

if [[ "${program}" == "entap" ]]; then
    cat $entaptsv | awk -F '\t' '{ print $1"\t"$25"\t"$28"\t"$29"\t"$32"\t"$2 }' >> simplifiedFunctional.tsv
else
    cat $entaptsv  | awk -F '\t' '{ print $1"\t"$2"\t"$5"\t"$8"\t"$9 }' >> simplifiedFunctional.tsv
fi
sed -i 's/ /_/g' simplifiedFunctional.tsv

sed -i '/^#/d' "$annotation"

python3 /seqprg/scripts/geneannotation/renameAnnots.py --annotation_file $annotation --functional_simplified simplifiedFunctional.tsv --output_file finalAnnotation.gtf

#rm the intermediate file
rm CorrectedAnnotation.gtf

agat_convert_sp_gxf2gxf.pl -g finalAnnotation.gtf -o finalAnnotation.gff3

sed -i '/AGAT\tfive_prime_UTR/d' finalAnnotation.gff3
sed -i '/AGAT\tthree_prime_UTR/d' finalAnnotation.gff3

rm finalAnnotation.gtf

# make the protein fasta file
agat_sp_extract_sequences.pl --clean_final_stop --gff finalAnnotation.gff3 -f genome.fa -p -o proteins_intermediate.fa

echo "Getting busco stats before running annotation checker"
# Get busco stats
conda activate BUSCO
busco -i proteins_intermediate.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate

rm proteins_intermediate.fa

mv buscoout/*/full_table.tsv .

rm -rf buscoout

mv finalAnnotation.gff3 semiFinalAnnotation.gff3

if [ -s "splign.gff3" ]; then
    echo "splign.gff3 exists and is populated."
    splign_command="--splign_file splign.gff3 "
else
    echo "splign.gff3 is missing or empty."
fi
if [ -s "reference_splign.gff3" ]; then
    echo "reference_splign.gff3 exists and is populated."
    reference_splign_command="--reference_splign_file reference_splign.gff3 "
else
    echo "reference_splign.gff3 is missing or empty."
fi
if [ -s "pasa_assembly.gff3" ]; then
    echo "pasa_assembly.gff3 exists and is populated."
    pasa_command="--pasa_file pasa_assembly.gff3 "
else
    echo "pasa_assembly.gff3 is missing or empty."
fi
if [ -s "reference_pasa_assembly.gff3" ]; then
    echo "reference_pasa_assembly.gff3 exists and is populated."
    reference_pasa_command="--reference_pasa_file reference_pasa_assembly.gff3 "
else
    echo "reference_pasa_assembly.gff3 is missing or empty."
fi
if [ -s "full_table.tsv" ]; then
    echo "full_table.tsv exists and is populated."
    busco_command="--busco_file full_table.tsv "
else
    echo "full_table.tsv is missing or empty."
fi
if [ -s "miniprot.gtf" ]; then
    echo "miniprot.gtf exists and is populated."
    miniprot_command="--protein_alignment_file miniprot.gtf "
else
    echo "miniprot.gtf is missing or empty."
fi
if [[ -z ${single_miniprot_coverage_min} ]]
then
    command_single_miniprot_coverage_min=""
else
    command_single_miniprot_coverage_min="--single_miniprot_coverage_min ${single_miniprot_coverage_min}"
fi
grep -E ';gene_biotype=(tRNA|pseudogene);' semiFinalAnnotation.gff3 > tRNAs.gff3
grep -v -E ';gene_biotype=(tRNA|pseudogene);' semiFinalAnnotation.gff3 > NotRNAs.gff3

rm semiFinalAnnotation.gff3

echo "RUNNING ANNOTATION CHECKER"
echo "python3 /seqprg/scripts/geneannotation/annotation_checker.py --annotation_file NotRNAs.gff3 ${splign_command}${reference_splign_command}${pasa_command}${reference_pasa_command}${busco_command}${miniprot_command}--multi_percent_same_exon_match ${multi_percent_same_exon_match} --multi_splign_coverage_min ${multi_splign_coverage_min} --multi_reference_splign_coverage_min ${multi_reference_splign_coverage_min} --multi_pasa_coverage_min ${multi_pasa_coverage_min} --multi_reference_pasa_coverage_min ${multi_reference_pasa_coverage_min} --multi_miniprot_coverage_min ${multi_miniprot_coverage_min} --multi_is_busco ${multi_is_busco} --multi_num_orthologs ${multi_num_orthologs} --single_percent_same_exon_match ${single_percent_same_exon_match} --single_splign_coverage_min ${single_splign_coverage_min} --single_reference_splign_coverage_min ${single_reference_splign_coverage_min} --single_pasa_coverage_min ${single_pasa_coverage_min} --single_reference_pasa_coverage_min ${single_reference_pasa_coverage_min} --single_is_busco ${single_is_busco} --single_num_orthologs ${single_num_orthologs} --target_ratio_single_multi {single_multi_ratio} --length_threshold ${length_threshold} ${command_single_miniprot_coverage_min}"
python3 /seqprg/scripts/geneannotation/annotation_checker.py --annotation_file NotRNAs.gff3 ${splign_command}${reference_splign_command}${pasa_command}${reference_pasa_command}${busco_command}${miniprot_command}--multi_percent_same_exon_match ${multi_percent_same_exon_match} --multi_splign_coverage_min ${multi_splign_coverage_min} --multi_reference_splign_coverage_min ${multi_reference_splign_coverage_min} --multi_pasa_coverage_min ${multi_pasa_coverage_min} --multi_reference_pasa_coverage_min ${multi_reference_pasa_coverage_min} --multi_miniprot_coverage_min ${multi_miniprot_coverage_min} --multi_is_busco ${multi_is_busco} --multi_num_orthologs ${multi_num_orthologs} --single_percent_same_exon_match ${single_percent_same_exon_match} --single_splign_coverage_min ${single_splign_coverage_min} --single_reference_splign_coverage_min ${single_reference_splign_coverage_min} --single_pasa_coverage_min ${single_pasa_coverage_min} --single_reference_pasa_coverage_min ${single_reference_pasa_coverage_min} --single_is_busco ${single_is_busco} --single_num_orthologs ${single_num_orthologs} --target_ratio_single_multi ${single_multi_ratio} --length_threshold ${length_threshold} ${command_single_miniprot_coverage_min}
#--splign_file splign.gff3 --reference_splign_file reference_splign.gff3 --pasa_file pasa_assembly.gff3 --reference_pasa_file reference_pasa_assembly.gff3 --busco_file full_table.tsv --protein_alignment_file miniprot.gtf

rm full_table.tsv

mv supported.gff3 finalAnnotation.gff3
mv removed.gff3 filtered_out_annotations.gff3
rm NotRNAs.gff3

cat tRNAs.gff3 >> finalAnnotation.gff3

rm tRNAs.gff3

agat_convert_sp_gxf2gxf.pl --gff finalAnnotation.gff3 -o finalAnnotationSorted.gff3

mv finalAnnotationSorted.gff3 finalAnnotation.gff3

agat_convert_sp_gff2gtf.pl --gff finalAnnotation.gff3 -o finalAnnotation.gtf

echo "adding comments to annotation"
currentTime=`date`
sed -i "1 i\#!Species: ${speciesName}" finalAnnotation.gtf
sed -i "1 i\#!Form Bio annotation workflow version 2.0" finalAnnotation.gtf
sed -i "1 i\#!Created with the Form Bio annotation workflow on: ${currentTime}" finalAnnotation.gtf
sed -i "1 i\##gtf-version 3" finalAnnotation.gtf

mv finalAnnotation.gtf final${speciesName}.gtf
mv finalAnnotation.gff3 final${speciesName}.gff3
sed -i 's/""""/""/g' final${speciesName}.gtf
sed -i 's/""""/""/g' final${speciesName}.gff3

sed -i '/\tfive_prime_UTR\t/d' final${speciesName}.gtf
sed -i '/\tthree_prime_UTR\t/d' final${speciesName}.gtf

sed -i '/\tfive_prime_UTR\t/d' final${speciesName}.gff3
sed -i '/\tthree_prime_UTR\t/d' final${speciesName}.gff3

# make the protein fasta file
echo "agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa -p -o proteins_${speciesName}.fa" >> parallel_agat_commands.txt

# make the cdna fasta file
echo "agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa --cdna -o cdna_${speciesName}.fa" >> parallel_agat_commands.txt

# make the mrna fasta file
echo "agat_sp_extract_sequences.pl --clean_final_stop --gff final${speciesName}.gff3 -f genome.fa --mrna -o mrna_${speciesName}.fa" >> parallel_agat_commands.txt

# convert gtf to gff to supply both in the final output
#agat_convert_sp_gxf2gxf.pl -g final${speciesName}.gtf -o final${speciesName}.gff3

# Get agat stats
echo "agat_sp_statistics.pl --gff final${speciesName}.gff3 --output final${speciesName}.AGAT.stats" >> parallel_agat_commands.txt

parallel < parallel_agat_commands.txt

sed -i 's/%25 /% /g' final${speciesName}.gff3

# Get busco stats
conda activate BUSCO
busco -i proteins_${speciesName}.fa -l ${lineage} -o buscoout -m protein -c ${threads}
conda deactivate

mv buscoout/short_summary.*.buscoout.txt .
cp short_summary.*.buscoout.txt busco.txt

# Make a file that only contains protein coding genes
sed '/gene_biotype "pseudogene"/d' final${speciesName}.gtf > final${speciesName}_nopseudo.gtf
sed '/gene_biotype "tRNA"/d' final${speciesName}_nopseudo.gtf > final${speciesName}_no_trna_or_pseudogenes.gtf

rm final${speciesName}_nopseudo.gtf

python3 /seqprg/scripts/geneannotation/finalReport.py --input_stats "final${speciesName}.AGAT.stats" --lineage "${lineage}" --scientific_name "${speciesName}"

conda activate biopython
python /seqprg/scripts/formatting/gff_to_genbank.py "final${speciesName}.gff3" "genome.fa"
conda deactivate
