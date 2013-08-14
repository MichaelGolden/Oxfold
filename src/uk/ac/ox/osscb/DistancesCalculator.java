package uk.ac.ox.osscb;

import uk.ac.ox.osscb.analysis.RNAFoldingTools;

public class DistancesCalculator {
	/**
	 * calculates distances between bases in alignment
	 * @author lepuslapis
	 * @param structure of alignment: if s[j]=k then j is paired with k; unpaired bases have s[i]=-1;
	 * @link {@link Constants#UnpairedBaseIdx} 
	 * @return distances between bases in alignment
	 */
	public static int[][] distCalc(int[] structure) {
		if (null==structure) {
			throw new IllegalArgumentException("Structure cannot be null.");
		}
		int length = structure.length;
		int[][] distances = new int[length][length];
		for (int b=1; b<length; b++) {
			for (int j=0; j<length-b; j++) {
				if ((structure[j]>j)&&(structure[j]<=j+b)) {
					distances[j][j+b]=1+distances[structure[j]][j+b];
				} else {
					//int tmp = distances[j+1][j+b];
					distances[j][j+b]=1+distances[j+1][j+b];
				}
			}
		}
		return distances;
	}
	
	public static int distCalc(int[] structure,int x, int y) {
		int[] distances = new int[y+1];
		for (int b=1; b<=y; b++) {
			int j = y - b;
			if ((structure[j]>j)&&(structure[j]<=j+b)) {
				distances[j]=1+distances[structure[j]];
			} else {
				distances[j]=1+distances[j+1];
			}
		}
		return distances[x];
	}
	
	public static int distCalc(int[] structure,int x, int y, int [][] distanceMatrix) {
		int[] distances = new int[y+1];
		for (int b=1; b<=y; b++) {
			int j = y - b;
			if ((structure[j]>j)&&(structure[j]<=j+b)) {
				distances[j]=1+distances[structure[j]];
				distanceMatrix[j][y] = distances[j];
			} else {
				distances[j]=1+distances[j+1];
				//distanceMatrix[j][y] = distances[j];
			}
		}
		return distances[x];
	}
	
	public static void main(String [] args)
	{
		//int [] pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString("((...))...(.......(...))(..(..))");
		//System.out.println(distCalc(pairedSites,0,pairedSites.length-1));
		//System.out.println(distCalc(pairedSites)[0][pairedSites.length-1]);
	}
}
