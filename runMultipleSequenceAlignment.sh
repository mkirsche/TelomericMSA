if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

usage() {
	echo ''
	echo 'Required arguments:'
	echo ''
	echo '  reference_genome   - fasta file containing the genome to align reads to'
	echo '  chromosome_name    - name of the chromosome (must match name in reference file)'
	echo '  telomere_length    - the length of telomere sequence which should be aligned '
	echo '  out_prefix         - prefix added to all output files (may include path as well)'
	echo '  reads_file         - fastq file containing sequencing reads'
	echo '  motif_start        - start position of any instance of the motif in the reference'
	echo '  motif_length       - length of the telomeric repeat motif'
	echo '  subtelomere_start  - start position of unique subtelomeric sequence in the reference'
	echo '  subtelomere_length - length of unique subtelomeric sequence to use'     
	echo ''
	echo 'Optional arguments:'
	echo ''
	echo '  readname_list      - file containing list of reads to target (one per line)'
	echo ''
}

while [ $# -gt 0 ]; do
  case "$1" in
    ref=* | reference_genome=*)
      REFERENCE_GENOME="${1#*=}"
      ;;
    chr=* | chromosome_name=*)
      CHRNAME="${1#*=}"
      ;;
    telomere_length=*)
      TELOMERE_LENGTH="${1#*=}"
      ;;
    out_prefix=*)
      OUTPREFIX="${1#*=}"
      ;;
    reads_file=*)
      READS_FILE="${1#*=}"
      ;;
    readname_list=*)
      READNAME_LIST="${1#*=}"
      ;;
    motif_start=*)
      MOTIF_START="${1#*=}"
      ;;
    motif_length=*)
      MOTIF_LENGTH="${1#*=}"
      ;;
    subtelomere_start=*)
      SUBTELOMERE_START="${1#*=}"
      ;;
    subtelomere_length=*)
      SUBTELOMERE_LENGTH="${1#*=}"
      ;;
    *)
      printf "Error: unknown option: $1\n"
      exit 1
  esac
  shift
done

missing_args=()

if [ -z $REFERENCE_GENOME ]
then
	missing_args+=("reference_genome")
fi

if [ -z $CHRNAME ]
then
	missing_args+=("chromosome_name")
fi

if [ -z $TELOMERE_LENGTH ]
then
	missing_args+=("telomere_length")
fi

if [ -z $OUTPREFIX ]
then
	missing_args+=("out_prefix")
fi

if [ -z $READS_FILE ]
then
	missing_args+=("reads_file")
fi

if [ -z $MOTIF_START ]
then
	missing_args+=("motif_start")
fi

if [ -z $MOTIF_LENGTH ]
then
	missing_args+=("motif_length")
fi

if [ -z $SUBTELOMERE_START ]
then
	missing_args+=("subtelomere_start")
fi

if [ -z $SUBTELOMERE_LENGTH ]
then
	missing_args+=("subtelomere_length")
fi

if [ ${#missing_args[@]} -gt 0 ]
then
  usage
  echo "Missing required arguments:"
  echo ''
  for i in "${missing_args[@]}"
  do
    echo '  '$i
  done
  echo ''
  exit 0
fi

samtools faidx $REFERENCE_GENOME
chrlength=`cat $REFERENCE_GENOME.fai | awk -v cn=$CHRNAME '{if($1 == cn){print $2}}'`
start=`expr $TELOMERE_LENGTH + 1`
#echo $start'-'$chrlength
CHROMOSOME_FILE=$OUTPREFIX'.chr.fa'
samtools faidx $REFERENCE_GENOME $CHRNAME:$start'-'$chrlength > $CHROMOSOME_FILE
PAF_FILE=$OUTPREFIX'.paf'
minimap2 -x map-ont $CHROMOSOME_FILE $READS_FILE > $PAF_FILE
SUBTELOMERE_END=`expr $SUBTELOMERE_START + $SUBTELOMERE_LENGTH - 1`
SUBMOTIF=`samtools faidx $REFERENCE_GENOME $CHRNAME:$SUBTELOMERE_START'-'$SUBTELOMERE_END | tail -n +2 | tr -d '\n' | tr ACGTacgt TGCATGCA | rev`

MOTIF_END=`expr $MOTIF_START + $MOTIF_LENGTH - 1`
MOTIF=`samtools faidx $REFERENCE_GENOME $CHRNAME:$MOTIF_START-$MOTIF_END | tail -n +2 | tr -d '\n' | tr ACGTacgt TGCATGCA | rev`

javac /home/mkirsche/eclipse-workspace/TelomereDetection/src/*.java
#echo $READS_FILE $PAF_FILE $OUTPREFIX $READNAME_LIST $MOTIF $SUBMOTIF
java -cp /home/mkirsche/eclipse-workspace/TelomereDetection/src MSA reads_file=$READS_FILE paf_file=$PAF_FILE out_prefix=$OUTPREFIX readlist_file=$READNAME_LIST motif=$MOTIF premotif=$SUBMOTIF

