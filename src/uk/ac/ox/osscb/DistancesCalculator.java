package uk.ac.ox.osscb;

public class DistancesCalculator {
	/**
	 * calculates distances between bases in alignment
	 * @author lepuslapis
	 * @param structure of alignment: if s[j]=k then j is paired with k; unpaired bases have s[i]=-1;
	 * @link {@link Constants#UnpairedBaseIdx} 
	 * @return distances between bases in alignment
	 */
	public int[][] distCalc(int[] structure) {
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
					int tmp = distances[j+1][j+b];
					distances[j][j+b]=1+distances[j+1][j+b];
				}
			}
		}
		return distances;
	}
}
