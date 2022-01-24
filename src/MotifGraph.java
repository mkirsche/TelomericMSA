/*
 * A graph representation of a telomeric motif repeat structure
 * 
 * The assumed structure is a unique pre-motif sequence followed by a particular
 * repeat sequence which is repeated an indefinite number of times.
 * 
 */
public class MotifGraph {
	
	// preMotif is the unique sequence, and is followed by the motif,
	// which is the repeated sequence
	String preMotif, motif;
	
	// The number of states in the graph is length(preMotif) + length(motif)
	int numStates;
	
	MotifGraph(String preMotif, String motif)
	{
		this.preMotif = preMotif;
		this.motif = motif;
		numStates = motif.length() + preMotif.length();
	}
	
	/*
	 * Find the optimal alignment path with a custom free length specific to this read
	 */
	String[] findOptimalPath(String s, int freeLength)
	{
		int tmp = Settings.freeLength;
		Settings.freeLength = freeLength;
		System.out.println(Settings.freeLength);
		String[] res = findOptimalPath(s);
		Settings.freeLength = tmp;
		return res;
	}
	
	/*
	 * Call dynamic programming function and then
	 * extract the actual alignment from traceback info
	 */
	String[] findOptimalPath(String s)
	{
		// Length of the read being aligned
		int m = s.length();
		
		// Used to store intermediate dynamic programming info to speed up algorithm
		memo = new double[m][numStates][3];
		
		// Keeps track of the optimal path during dynamic programming
		trace = new int[m][numStates][3];
		
		// Call the dynamic programming alignment on every sub-problem to fill optimal path info
		for(int i = m; i>=0; i--)
		{
			for(int j = 0; j<preMotif.length(); j++)
			{
				for(int k = 0; k<3; k++)
				{
					bestScore(s, i, j+motif.length(), k);
				}
			}
			for(int j = 0; j<motif.length(); j++)
			{
				for(int k = 0; k<3; k++)
				{
					bestScore(s, i, j, k);
				}
			}
		}
		
		// The index in the read sequence
		int index = 0;
		
		// The node in the graph (initially the start of the pre-motif
		int state = motif.length();
		
		// Whether the last action was a match/mismatch (0), insertion (+), or deletion (-)
		int indel = 0;
		
		// The two aligned strings
		String alignedRef = "";
		String alignedRead = "";
		while(index < s.length())
		{
			int val = trace[index][state][indel];
			if(val <= -2*motif.length())
			{
				// Deletion before entering the motif sequence
				int numDeleted = (Math.abs(val) - 2*motif.length())%preMotif.length();
				for(int i = 0; i<numDeleted; i++)
				{
					alignedRef += preMotif.charAt(state - motif.length() + i);
					alignedRead += "-";
				}
				state += numDeleted;
				if(state == motif.length() || state == numStates) state = 0;
				if(Math.abs(val) >= preMotif.length())
				{
					// Insertion after deletion
					alignedRef+= "-";
					alignedRead += s.charAt(index);
					index++;
					indel = 1;
				}
				else
				{
					// Match/mismatch after deletion
					alignedRef += preMotif.charAt(state - motif.length());
					alignedRead += s.charAt(index);
					state++;
					if(state == motif.length() || state == numStates) state = 0;
					index++;
					indel = 0;
				}
				continue;
			}
			if(val == 0)
			{
				// Match/mismatch
				alignedRef += state >= motif.length() ? preMotif.charAt(state - motif.length()) : motif.charAt(state);
				alignedRead += s.charAt(index);
				index++;
				state++;
				if(state == motif.length() || state == numStates) state = 0;
			}
			else if(val > 0)
			{
				// Insertion
				alignedRef += "-";
				alignedRead += s.charAt(index);
				index++;
				indel = 1;
			}
			else if(val < 0)
			{
				// Deletion after entering motif sequence
				val *= -1;
				int numDeleted = val%motif.length();
				for(int i = 0; i<numDeleted; i++)
				{
					alignedRef += motif.charAt((state + i)%motif.length());
					alignedRead += "-";
				}
				state = (state + numDeleted)%motif.length();
				if(val >= motif.length())
				{
					// Insertion after deletion
					alignedRef += "-";
					alignedRead += s.charAt(index);
					index++;
					indel = 1;
				}
				else
				{
					// Match/mismatch after deletion
					alignedRef += motif.charAt(state);
					alignedRead += s.charAt(index);
					state++;
					if(state == motif.length() || state == numStates) state = 0;
					index++;
					indel = 0;
				}
			}
		}
		return new String[] {alignedRef, alignedRead};
	}
	
	/*
	 * indel is 0 if not currently inserting/deleting, 1 if inserting, 2 if deleting
	 */
	double[][][] memo;
	int[][][] trace;

	double bestScore(String s, int index, int state, int indel)
	{
		if(index == s.length())
		{
			return 1;
		}
		
		if(state == numStates)
		{
			return bestScore(s, index, 0, indel);
		}
		
		if(memo[index][state][indel] != 0)
		{
			return memo[index][state][indel];
		}
		
		if(state >= motif.length())
		{
			// Still in pre-motif so no cycling happens
			/*
			 * Cases:
			 *   1. Advance state and position: match or mismatch depending on character equality
			 *   2. Advance only state (deletion) - 
			 *      try all possible new states but also follow with case 1 to avoid cycles
			 *   3. Advance only position (insertion)
			 */
			double[] probs = new double[4];
			if(indel == 0)
			{
				probs[1] = Settings.probSubstitution;
				probs[2] = Settings.probInsertion;
				probs[3] = Settings.probDeletion;
			}
			else if(indel == 1)
			{
				probs[1] = Settings.probSubstitution;
				probs[2] = Settings.probInsertionContinue;
				probs[3] = Settings.probDeletion;
			}
			else if(indel == 2)
			{
				probs[1] = Settings.probSubstitution;
				probs[2] = Settings.probInsertion;
				probs[3] = Settings.probDeletionContinue;
			}
			probs[0] = 1 - probs[1] - probs[2] - probs[3];
			
			double res = -1e9;
			
			// Match or mismatch
			if(state != numStates-1)
			{
				double case1 = Math.log(preMotif.charAt(state - motif.length()) == s.charAt(index) ? probs[0] : probs[1])
					+ bestScore(s, index+1, state + 1, 0);
				res = Math.max(res, case1);
				trace[index][state][indel] = 0;
			}
			
			// Some amount of deletion
			double runningProb = Math.log(Settings.probDeletion);
			for(int i = 0; i<preMotif.length(); i++)
			{
				if(state + i + 2 <= numStates)
				{
					double case2Match = runningProb 
							+ bestScore(s, index + 1, (state + i + 2), 0) 
							+ Math.log(preMotif.charAt(i+state+1-motif.length()) == s.charAt(index) ? probs[0] : probs[1]);
					if(index+ Settings.freeLength < s.length() && case2Match > res)
					{
						res = case2Match;
						trace[index][state][indel] = -(i + 1) - 2*motif.length();
					}
				}
				if(state + i + 1 <= numStates)
				{
					double case2Insertion = runningProb + Math.log(Settings.probInsertion) + bestScore(s, index + 1, (state + i + 1), 1);
					if(index+ Settings.freeLength < s.length() && case2Insertion > res)
					{
						res = case2Insertion;
						trace[index][state][indel] = -(i + 1) - preMotif.length() - 2*motif.length();
					}
				}
				else
				{
					break;
				}
				
				runningProb += Math.log(Settings.probDeletionContinue);
			}
			double case3 = Math.log(probs[2]) + bestScore(s, index + 1, state, 1);
			if(index+ Settings.freeLength < s.length() && case3 > res)
			{
				res = case3;
				trace[index][state][indel] = 1;
			}
			return memo[index][state][indel] = res;
		}
				
		/*
		 * Cases:
		 *   1. Advance state and position: match or mismatch depending on character equality
		 *   2. Advance only state (deletion) - 
		 *      try all possible new states but also follow with case 1 to avoid cycles
		 *   3. Advance only position (insertion)
		 */
		// Probabilities of match/substitution/insertion/deletion
		double[] probs = new double[4];
		if(indel == 0)
		{
			// Last update was match/mismatch
			probs[1] = Settings.probSubstitution;
			probs[2] = Settings.probInsertion;
			probs[3] = Settings.probDeletion;
		}
		else if(indel == 1)
		{
			// Last update was insertion
			probs[1] = Settings.probSubstitution;
			probs[2] = Settings.probInsertionContinue;
			probs[3] = Settings.probDeletion;
		}
		else if(indel == 2)
		{
			//Last update was deletion
			probs[1] = Settings.probSubstitution;
			probs[2] = Settings.probInsertion;
			probs[3] = Settings.probDeletionContinue;
		}
		probs[0] = 1 - probs[1] - probs[2] - probs[3];
		
		// Case 1: Match or mismatch
		double case1 = Math.log(motif.charAt(state) == s.charAt(index) ? probs[0] : probs[1])
				+ bestScore(s, index+1, (state + 1) % motif.length(), 0);
		double res = case1;
		
		// Trace value of 0 - match or mismatch
		trace[index][state][indel] = 0;
		
		// Case 2: Delete one or more bases
		// Must include a non-deletion at the end to avoid infinite loop because of cyclic motif
		
		// Keep track of running probability for multi-base deletions
		double runningProb = Math.log(Settings.probDeletion);
		
		// Try deleting anywhere from 1 to len(motif) bases
		for(int i = 0; i<motif.length(); i++)
		{
			// Deletion followed by match/mismatch
			double case2Match = runningProb 
					+ bestScore(s, index + 1, (state + i + 2)%motif.length(), 0) 
					+ Math.log(motif.charAt((i+state+1)%motif.length()) == s.charAt(index) ? probs[0] : probs[1]);
			if(index+ Settings.freeLength < s.length() && case2Match > res)
			{
				res = case2Match;
				// Negative trace value < len(motif) is deletion then match
				trace[index][state][indel] = -(i + 1);
			}
			
			// Deletion followed by inserted base
			double case2Insertion = runningProb + Math.log(Settings.probInsertion) + bestScore(s, index + 1, (state + i + 1)%motif.length(), 1);
			if(index + Settings.freeLength < s.length() && case2Insertion > res)
			{
				res = case2Insertion;
				// Negative trace value in [len(motif),2*len(motif)) is deletion then insertion
				trace[index][state][indel] = -(i + 1) - motif.length();
			}
			
			// Update running probability
			runningProb += Math.log(Settings.probDeletionContinue);
		}
		
		// Case 3: Insert a base
		double case3 = Math.log(probs[2]) + bestScore(s, index + 1, (state)%motif.length(), 1);
		if(index+ Settings.freeLength < s.length() && case3 > res)
		{
			res = case3;
			// Positive trace value indicates insertion
			trace[index][state][indel] = 1;
		}
		return memo[index][state][indel] = res;
	}
}
