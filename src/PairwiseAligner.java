import java.util.ArrayList;

public class PairwiseAligner {

	int mismatchPenalty;
	int matchScore;
	int gapStartPenalty;
	int gapContinuePenalty;
	
	PairwiseAligner(int mm, int ms, int ss, int cc)
	{
		mismatchPenalty = mm;
		matchScore = ms;
		gapStartPenalty = ss;
		gapContinuePenalty = cc;
	}
	
	Alignment findBest(String s, String t)
	{
		int n = s.length();
		int m = t.length();
		// Third dimension is gap state - 0 means no gap, 1 means gap in s, 2 means gap in t
		int[][][] scores = new int[n+1][m+1][3];
		int[][][] next = new int[n+1][m+1][3];
		for(int i = n; i>=0; i--)
		{
			for(int j = m; j>=0; j--)
			{
				for(int k = 0; k<3; k++)
				{
					if(i == n && j == m)
					{
						scores[i][j][k] = 0;
						continue;
					}
					scores[i][j][k] = (int)-1e9;
					if(i < n)
					{
						// Gap in t
						int newScore = scores[i+1][j][2] + (k == 2 ? gapContinuePenalty : gapStartPenalty);
						if(newScore > scores[i][j][k])
						{
							scores[i][j][k] = newScore;
							next[i][j][k] = 2;
						}
					}
					if(j < m)
					{
						// Gap in 2
						int newScore = scores[i][j+1][1] + (k == 1 ? gapContinuePenalty : gapStartPenalty);
						if(newScore > scores[i][j][k])
						{
							scores[i][j][k] = newScore;
							next[i][j][k] = 1;
						}
					}
					if(i < n && j < m)
					{
						int newScore = scores[i+1][j+1][k] + (s.charAt(i) == t.charAt(j) ? matchScore : mismatchPenalty);
						if(newScore > scores[i][j][k])
						{
							scores[i][j][k] = newScore;
							next[i][j][k] = 0;
						}
					}
				}
			}
		}
		return new Alignment(scores, next);
	}
	
	static class Alignment
	{
		int[] vector;
		int score;
		Alignment(int[][][] scores, int[][][] next)
		{
			ArrayList<Integer> steps = new ArrayList<Integer>();
			score = scores[0][0][0];
			int i = 0, j = 0, k = 0;
			while(i < scores.length - 1 || j < scores[0].length - 1)
			{
				steps.add(next[i][j][k]);
				if(next[i][j][k] == 0)
				{
					i++;
					j++;
					k = 0;
				}
				else if(next[i][j][k] == 1)
				{
					j++;
					k = 1;
				}
				else if(next[i][j][k] == 2)
				{
					i++;
					k = 2;
				}
			}
			vector = new int[steps.size()];
			for(i = 0; i<vector.length; i++)
			{
				vector[i] = steps.get(i);
			}
		}
	}
}
