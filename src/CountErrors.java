import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.Scanner;

public class CountErrors {
public static void main(String[] args) throws Exception
{
	/*if(args.length != 2)
	{
		System.out.println("Usage: java CountErrors [input aln file] [output tsv file]");
		System.exit(0);
	}*/
	String fn = args[0];//"/home/mkirsche/eclipse-workspace/TelomereDetection/bigtest_msa.aln";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	HashMap<String, String> nameToAlnSeq = new HashMap<String, String>();
	HashMap<String, ArrayDeque<Integer>> errors = new HashMap<String, ArrayDeque<Integer>>();
	PrintWriter out = new PrintWriter(new File(args[1]));//"bigtest_counts.txt"));
	while(input.hasNext())
	{
		String line = input.nextLine();
		
		// Empty line here means the end of an alignment block, so process the current block
		if(line.length() == 0)
		{
			if(nameToAlnSeq.size() > 0)
			{
				updateErrors(nameToAlnSeq, errors, "template");
				nameToAlnSeq = new HashMap<String, String>();
			}
		}
		else
		{
			//Remove any occurrences of multiple spaces for easier tokenizing
			while(line.contains("  "))
			{
				line = line.replaceAll("  ", " ");
			}
			
			// The line should have a read name and an alignment String separated by a space
			String[] tokens = line.trim().split(" ");
			if(tokens.length != 2)
			{
				continue;
			}
			
			// Add this line to the name->alignment map
			nameToAlnSeq.put(tokens[0], tokens[1]);
		}
	}
	
	// If there are leftover alignments from the last block, process those
	if(nameToAlnSeq.size() > 0)
	{
		updateErrors(nameToAlnSeq, errors, "template");
		nameToAlnSeq = new HashMap<String, String>();
	}
	
	// Write results to tsv - X is distance from read end and Y is number of errors
	out.println("READNAME\tX\tY");
	for(String s : errors.keySet())
	{
		ArrayDeque<Integer> list = errors.get(s);
		while(!list.isEmpty())
		{
			out.println(s+"\t"+list.size()+"\t"+list.pollFirst());
		}
	}
	out.close();
}

/*
 * Updates the error counts of all sequences based on a block of the MSA file
 * 
 * nameToAlnSeq
 */
static void updateErrors(HashMap<String, String> nameToAlnSeq, HashMap<String, ArrayDeque<Integer>> errors, String refName)
{
	// The number of aligned positions in this block
	int length = nameToAlnSeq.get(refName).length();
	
	for(int i = 0; i<length; i++)
	{
		// The current character of the reference - used to detect mismatches and ref position updates
		char refAlignmentChar = nameToAlnSeq.get(refName).charAt(i);
		
		for(String s : nameToAlnSeq.keySet())
		{
			if(s.equals(refName))
			{
				continue;
			}
			
			String currentAlnSeq = nameToAlnSeq.get(s);
			char currentAlnChar = currentAlnSeq.charAt(i);
			if(currentAlnChar == 'x')
			{
				continue;
			}
			boolean mismatch = false;
			if(currentAlnChar != 'x' && currentAlnChar != refAlignmentChar)
			{
				mismatch = true;
			}
			
			if(!errors.containsKey(s))
			{
				errors.put(s, new ArrayDeque<Integer>());
			}
			ArrayDeque<Integer> errorCountList = errors.get(s);
			if(errorCountList.size() == 0)
			{
				errorCountList.add(0);
			}
			
			if(refAlignmentChar == '-')
			{
				errorCountList.add(errorCountList.pollLast() + (mismatch ? 1 : 0));
			}
			else
			{
				errorCountList.add(errorCountList.peekLast() + (mismatch ? 1 : 0));
			}
		}
		
		if(refAlignmentChar != '-')
		{
			
		}
	}
}
}
