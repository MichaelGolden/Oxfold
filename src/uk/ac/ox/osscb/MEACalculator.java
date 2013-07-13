package uk.ac.ox.osscb;

import java.math.BigDecimal;
import java.util.Arrays;

import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;

public class MEACalculator {
	
		public int[] calculate(PosteriorProbabilities posteriorProbabilities) {
			System.err.println("WARNING: MEACalculator may be broken, sometimes predicts impossible base-pairs");
			BigDecimal[][] pairedProbs = posteriorProbabilities.getPairedProbs();
			BigDecimal[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
			int length = unpairedProbs.length;
						
			BigDecimal[][] eMatrix = new BigDecimal[length][length];
			int[][] dMatrix = new int[length][length]; //used for decoding
			
			for (int j = 0; j<length; j++) {
				BigDecimal tmp = unpairedProbs[j].multiply(BigDecimal.valueOf(0.5));
				eMatrix[j][j] = tmp; dMatrix[j][j] = j;
			}
			for (int b = 1; b<length; b++) {
				for (int j = 0; j<length-b; j++) {
					eMatrix[j][j+b] = eMatrix[j][j].add(eMatrix[j+1][j+b]);
					dMatrix[j][j+b] = j;
				}
			}
						
			/*
			 * do the MAP recursion
			 */
			for (int b = 3; b<length; b++) {
				for (int j = 0; j<length-b; j++) {
					for (int l = j+2; l<j+b; l++) {
						BigDecimal tmp = pairedProbs[j][l].add(eMatrix[j+1][l-1]).add(eMatrix[l+1][j+b]).subtract(Constants.MEApenalty);;
						if (tmp.compareTo(eMatrix[j][j+b])>0) {
							eMatrix[j][j+b] = tmp; dMatrix[j][j+b] = l;
						}
					}
					BigDecimal tmp = pairedProbs[j][j+b].add(eMatrix[j+1][j+b-1]);
					if (tmp.compareTo(eMatrix[j][j+b])>0) {
						eMatrix[j][j+b] = tmp; dMatrix[j][j+b] = j+b;
					}
				}
			}
			
			/*
			 * do the decoding (go through a stack)
			 */
			int[] structure = new int[length];
			Arrays.fill(structure, Constants.UnpairedBaseIdx);
			decode(structure,dMatrix);		
 			
			return structure;
		}
		
		private void decode(int[] structure, int[][] dMatrix) {
			int length = structure.length;
			int[] left = new int[length]; int[] right = new int[length];
			left[0] = 0; right[0] = length-1;
			int last = 0;
			while (last>=0) {
				int l = left[last]; int r = right[last];
				if (l >= r) {
					structure[r] = -1;
					last--;
				} else {
					int d = dMatrix[l][r];
					if (d == l) {
						structure[l] = -1;
						left[last] = l+1;
					} else if (d == r) {
						structure[l] = r; structure[r] = l;
						left[last] = l+1;
						right[last] = r-1;
					} else {
						structure[l] = d; structure[d] = l;
						left[last] = l+1; right[last] = d-1;
						left[last+1] = d+1; right[last+1] = r;
						last++;	
					}
				}
			}
		}
		
}
