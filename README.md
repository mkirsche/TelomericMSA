# Telomeric Multiple Sequence Alignment
Script for performing anchored free-end multiple sequence alignment on telomeric sequences for detected elongation.


## Software Dependencies

[samtools](https://anaconda.org/bioconda/samtools) and [minimap2](https://anaconda.org/bioconda/minimap2) are both required to run this software. It is assumed they are on the user's path, but if not, the calls to these tools in runMultipleSequenceAlignment.sh can be updated to point to a specific installation.

In addition, the source code for this software is in Java, so a Java Development Kit (jdk) is required, such as [openjdk](https://anaconda.org/conda-forge/openjdk) which is available through conda.


## Usage

```
./runMultipleSequenceAlignment.sh <args>

Required arguments:

  reference_genome   - fasta file containing the genome to align reads to
  chromosome_name    - name of the chromosome (must match name in reference file)
  telomere_length    - the length of telomere sequence which should be aligned 
  out_prefix         - prefix added to all output files (may include path as well)
  reads_file         - fastq file containing sequencing reads
  motif_start        - start position of any instance of the motif in the reference
  motif_length       - length of the telomeric repeat motif
  subtelomere_start  - start position of unique subtelomeric sequence in the reference
  subtelomere_length - length of unique subtelomeric sequence to use

Optional arguments:

  readname_list      - file containing list of reads to target (one per line)
  free_length        - number of bases at the end of each read which is aligned freely (default 50)
  global_free_length - number of bases from the end of longest read which is aligned freely (default -1, meaning no global free alignment)

```

Arguments should be passed using the format `<argument_name>=value`. Running this script with no parameters will give a full usage menu.  An example of running it can be found in `run_sample.sh`.



