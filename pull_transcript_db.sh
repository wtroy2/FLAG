# This script can be used to pull a blast transcript database
## Make sure that you already have the singularity or docker images to do this
usage() {
  echo "-h Help documentation for pull_transcript_db.sh"
  echo "-c  --containerization software (docker or singularity)"
  echo "Example: bash pull_transcript_db.sh -c docker"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :c:h opt
do
    case $opt in
        c) containerization=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z ${containerization} ]]
then
    containerization="docker"
fi

cd databases
if [[ "$containerization" == "docker" ]]; then
    echo "Using docker to pull refseq_select_rna"
    echo "docker run -it -v $(pwd):/data ghcr.io/formbio/flag_ncbiclibraries:latest perl /opt/ncbi-blast-2.13.0+/bin/update_blastdb.pl refseq_select_rna"
    docker run -it -v $(pwd):/data ghcr.io/formbio/flag_ncbiclibraries:latest perl /opt/ncbi-blast-2.13.0+/bin/update_blastdb.pl refseq_select_rna
elif [[ "$containerization" == "singularity" ]]; then
    echo "Using singularity to pull refseq_select_rna."
    echo "Please make sure your singularity env is active if you are having trouble or just use the command below in your terminal."
    echo "singularity exec --bind $(pwd):/data ../containers/ncbiclibraries/flag_ncbiclibraries.image perl /opt/ncbi-blast-2.13.0+/bin/update_blastdb.pl refseq_select_rna"
    singularity exec --bind $(pwd):/data ../containers/ncbiclibraries/flag_ncbiclibraries.image perl /opt/ncbi-blast-2.13.0+/bin/update_blastdb.pl refseq_select_rna
else
    echo "Please choose to use docker or singularity"
fi
cd ..
