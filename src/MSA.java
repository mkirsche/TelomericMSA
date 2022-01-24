import java.util.*;
import java.io.*;
public class MSA {
public static void main(String[] args) throws Exception
{
	if(args.length == 0)
	{
		args = new String[] {
				"reads_file=/home/mkirsche/telomere/single_1L_telomere_reads.fastq",
				"paf_file=/home/mkirsche/telomere/alignments.paf",
				"out_prefix=/home/mkirsche/telomere/test",
				"motif=TGGTGTGTGGGTG",
				"premotif=CGAAGGCTTTAATTTGCCGGAGTACTGTCCTCCGAGCGGTGAACTGTCCTCCGACTCGTGCGGAGTACTGTCCTCCGATCGGTGAACTGTCCTCCGCGAATTCCGGAGTACTAATCTCCG",
				"readlist_file=/home/mkirsche/telomere/msa_gappenalty5.aln"
		};
	}
	Settings.parseArgs(args);
	HashSet<String> readNames = extractReadNames();

	Scanner input = new Scanner(new FileInputStream(new File(Settings.readsFn)));
	
	HashMap<String, String> readNameToSeq = new HashMap<String, String>();
	
	String name = "";
	StringBuilder seq = new StringBuilder("");
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.charAt(0) == '@')
		{
			if(name.length() > 0)
			{
				String key = getKey(readNames, name);
				if(key.length() > 0)
				{
					readNameToSeq.put(key+"_"+seq.toString().length(), seq.toString());
				}
			}
			seq = new StringBuilder("");
			name = line.substring(1).split(" ")[0];
			continue;
		}
		else if(line.equals("+"))
		{
			input.nextLine();
			continue;
		}
		
		seq.append(line);
	}
	if(name.length() > 0)
	{
		String key = getKey(readNames, name);
		if(key.length() > 0)
		{
			readNameToSeq.put(key+"_"+seq.toString().length(), seq.toString());
		}
	}
	
	input = new Scanner(new FileInputStream(new File(Settings.pafFile)));
	HashMap<String, String> readNameToTelomere = new HashMap<String, String>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		String readName = tokens[0];
		String key = getKey(readNames, readName);
		if(key.length() > 0)
		{
			int readEnd = Integer.parseInt(tokens[3]);
			char strand = tokens[4].charAt(0);
			String totalSeq = readNameToSeq.get(key+"_"+Integer.parseInt(tokens[1]));
			String telomereSeq = "";
			if(strand == '-')
			{
				if(readEnd < 120) continue;
				telomereSeq = revComp(totalSeq.substring(readEnd - 120));
			}
			else
			{
				continue;
			}
			if(telomereSeq.length() > 2000)
			{
				continue;
			}
			readNameToTelomere.put(key, telomereSeq);
		}
	}
	
	int maxLength = 0;
	for(String s : readNameToTelomere.keySet())
	{
		String telomereSeq = readNameToTelomere.get(s);
		if(maxLength == 0)
		{
			maxLength = telomereSeq.length();
		}
		maxLength = Math.max(maxLength, telomereSeq.length());
	}
	
	MotifGraph mg = new MotifGraph(Settings.preMotif, Settings.motif);
	
	String[] names = new String[readNameToTelomere.size()];
	String[][] alignments = new String[readNameToTelomere.size()][];
	int idx = 0;
	int readsAligned = 0;
	for(String s : readNameToTelomere.keySet())
	{
		if(readsAligned > 0 && readsAligned%10 == 0)
		{
			System.err.println("Aligned " + readsAligned + " out of " + readNameToTelomere.size());
		}
		names[idx] = s;
		String telomereSeq = readNameToTelomere.get(s);
		telomereSeq = revComp(telomereSeq);
		if(Settings.globalFreeLength == -1)
		{
			alignments[idx] = mg.findOptimalPath(telomereSeq);
		}
		else
		{
			alignments[idx] = mg.findOptimalPath(telomereSeq, Math.max(Settings.freeLength, Settings.globalFreeLength + telomereSeq.length() - maxLength));
		}
		idx++;
		readsAligned++;
	}
	
	printAlignments(alignments, names, false, Settings.outPrefix + "_msa_full.aln");
	printAlignments(alignments, names, true, Settings.outPrefix + "_msa.aln");
}

/*
 * Extract the names of reads that need to be considered
 */
static HashSet<String> extractReadNames() throws Exception
{
	HashSet<String> readNames = new HashSet<String>();
	if(Settings.readsListFile.length() == 0)
	{
		return readNames;
	}
	Scanner input = new Scanner(new FileInputStream(new File(Settings.readsListFile)));
	//input.nextLine();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.charAt(0) == ' ' || line.charAt(0) == '\t') continue;
		String readName = line.split(" ")[0];
		if(readNames.contains(readName)) break;
		readNames.add(readName);
	}
	return readNames;
}

/*
 * 
 */
static void printAlignments(String[][] alignments, String[] readNames, boolean filterSeqError, String outfile) throws Exception
{
	String[] res = combine(alignments, filterSeqError);
	
	PrintWriter out = new PrintWriter(new File(outfile));
	out.println("Multiple sequence alignment");
	out.println();
	out.println();
	int length = res[0].length();
	int maxName = 0;
	for(String s : readNames) maxName = Math.max(maxName, s.length());
	String template = "template";
	for(int i = 0; i<maxName - 8; i++) template += " ";
	for(int i = 0; i<length; i+=50)
	{
		for(int j = 0; j<res.length; j++)
		{
			out.println((j < res.length - 1 ? readNames[j] : template) +"      " + res[j].substring(i, Math.min(res[j].length(), i+50)));
		}
		out.println();
		out.println();
	}
	out.close();
}

/*
 * Get the name of the read from the set of important reads that
 * either matches or is a prefix of name
 */
static String getKey(HashSet<String> readNames, String str)
{
	if(readNames.size() == 0)
	{
		return str;
	}
	String res = "";
	for(String s : readNames)
	{
		if(s.equals(str))
		{
			return s;
		}
		if(str.startsWith(s))
		{
			return s;
		}
	}
	return res;
}

static String[] combine(String[][] alignments, boolean filterErrors)
{
	int n = alignments.length;
	StringBuilder[] preRes = new StringBuilder[n+1];
	for(int i = 0; i<n+1; i++)
	{
		preRes[i] = new StringBuilder("");
	}
	int[] ptrs = new int[n];
	while(true)
	{
		boolean done = true;
		boolean blank = false;
		boolean nearEnd = filterErrors ? false : true;
		for(int i = 0; i<n; i++)
		{
			done &= ptrs[i] == alignments[i][0].length();
			if(ptrs[i] != alignments[i][0].length())
			{
				blank |= alignments[i][0].charAt(ptrs[i]) == '-';
				nearEnd |= ptrs[i] + 100 >= alignments[i][0].length();
			}
			else
			{
				nearEnd = true;
			}
		}
		if(done)
		{
			break;
		}
		
		if(blank)
		{
			if(nearEnd) preRes[n].append('-');
			for(int i = 0; i<n; i++)
			{
				if(ptrs[i] != alignments[i][0].length())
				{
					if(alignments[i][0].charAt(ptrs[i]) == '-')
					{
						if(nearEnd) preRes[i].append(alignments[i][1].charAt(ptrs[i]));
						ptrs[i]++;
					}
					else
					{
						if(nearEnd) preRes[i].append('-');
					}
				}
				else
				{
					if(nearEnd) preRes[i].append('x');
				}
			}
		}
		else
		{
			char c = '-';
			for(int i = 0; i<n; i++)
			{
				if(ptrs[i] != alignments[i][0].length())
				{
					c = alignments[i][0].charAt(ptrs[i]);
					preRes[i].append(alignments[i][1].charAt(ptrs[i]));
					ptrs[i]++;
				}
				else
				{
					preRes[i].append('x');
				}
			}
			preRes[n].append(c);
		}
	}
	String[] res = new String[n+1];
	for(int i = 0; i<n+1; i++)
	{
		res[i] = preRes[i].toString();
	}
	return res;
	
}

static String revComp(String s)
{
	char[] res = new char[s.length()];
	for(int i = 0; i<s.length(); i++)
	{
		char c = s.charAt(s.length() - 1 - i);
		res[i] = c;
		if(c == 'A') res[i] = 'T';
		else if(c == 'C') res[i] = 'G';
		else if(c == 'G') res[i] = 'C';
		else if(c == 'T') res[i] = 'A';
	}
	return new String(res);
}

static class Read
{
	String name, seq;
	Read(String name, String seq)
	{
		this.name = name;
		this.seq = seq;
	}
}
}
