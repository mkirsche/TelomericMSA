if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

REFERENCE_GENOME=$1
CHRNAME=$2
TELOMERE_LENGTH=$3
OUTPREFIX=$4
READS_FILE=$5
READNAME_LIST=$6
MOTIF_START=$7
MOTIF_LENGTH=$8
SUBTELOMERE_START=$9
SUBTELOMERE_LENGTH=${10}

#CHRNAME=chrI
#READS_FILE=
#READNAME_LIST=/home/mkirsche/telomere/msa_gappenalty5.aln
#MOTIF=TGGTGTGTGGGTG
#MOTIF_START=1027
#MOTIF_LENGTH=13
#SUBTELOMERE_START=1066
#SUBTELOMERE_LENGTH=120

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

