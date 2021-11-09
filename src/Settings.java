
public class Settings {
	static String preMotif = "";
	static String motif = "";
	static String readsListFile = "";
	static String readsFn = "";
	static String outPrefix = "";
	static String pafFile = "";
	static int freeLength = 50;
	static int subtelomereLength = 120;
	
	static double probInsertion = .02;
	static double probInsertionContinue = .1;
	static double probDeletionContinue = .2;
	static double probDeletion = .05;
	static double probSubstitution = .01;
	
	/*
	 * Print the usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Free-end anchored telomeric multiple-sequence aligner");
		System.out.println("Usage: java -cp src MSA [args]");
		System.out.println("  Example: java -cp src MSA file_list=filelist.txt out_file=out.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  reads_file (String) - a FASTQ file containing the reads to align");
		System.out.println("  paf_file   (String) - a PAF file containing alignments of the reads to the chromosome");
		System.out.println("  out_prefix (String) - the prefix for the .aln files which are produced");
		System.out.println("  motif      (String) - the telomeric motif which is repeated to form telomeric sequences");
		System.out.println("  premotif   (String) - unique sub-telomeric sequence used to anchor the alignments");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  free_length           (int)    [50] - the length of sequence at ends of reads which is not able to contain gaps");
		System.out.println("  subtelomere_length    (int)    [50] - the length of unique subtelomeric sequence to keep from each read");
		System.out.println("  readlist_file  (String)      - a file contaning list of read names to consider, one per line");
		System.out.println();
	}
	
	/*
	 * Method for parsing all command line arguments and making sure the required
	 * arguments are given values
	 */
	static void parseArgs(String[] args) throws Exception
	{
		if(args.length == 1 && (args[0].equalsIgnoreCase("--version") || args[0].equalsIgnoreCase("-v")))
		{
			System.out.println("Free-end anchored telomeric MSA version 1.0.0");
			System.exit(0);
		}
		if(args.length == 0)
		{
			usage();
			System.exit(0);
		}
		
		for(int i = 0; i<args.length; i++)
		{
			if(args[i].indexOf('=') == -1)
			{
				
			}
			else
			{
				int equalIdx = args[i].indexOf('=');
				String key = args[i].substring(0, equalIdx);
				while(key.length() > 0 && key.charAt(0) == '-')
				{
					key = key.substring(1);
				}
				String val = args[i].substring(1 + equalIdx);
				
				switch(key)
				{
					case("reads_file"):
						readsFn = val;
						break;
					case("paf_file"):
						pafFile = val;
						break;
					case("readlist_file"):
						readsListFile = val;
						break;
					case("out_prefix"):
						outPrefix = val;
						break;
					case("motif"):
						motif = val;
						break;
					case("premotif"):
						preMotif = val;
						break;
					case("free_length"):
						freeLength = Integer.parseInt(val);
						break;
					case("subtelomere_length"):
						subtelomereLength = Integer.parseInt(val);
						break;
					default:
						break;
				}
			}
		}
		
		if(readsFn.length() == 0 || readsListFile.length() == 0 || outPrefix.length() == 0 || 
				motif.length() == 0 || preMotif.length() == 0 || pafFile.length() == 0)
		{
			usage();
			System.exit(0);
		}
			
	}
	
}
