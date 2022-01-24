if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

$BINDIR/runMultipleSequenceAlignment.sh reference_genome=$BINDIR/reference_genome.fa chromosome_name=chrI telomere_length=1000 out_prefix=$BINDIR/bigtest reads_file=$BINDIR/single_1L_telomere_reads.fastq motif_start=1027 motif_length=13 subtelomere_start=1066 subtelomere_length=120 free_length=0 global_free_length=800
