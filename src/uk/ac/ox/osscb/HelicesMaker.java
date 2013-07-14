package uk.ac.ox.osscb;



import uk.ac.ox.osscb.inoutside.Helix;

/**
 * 
 * @author Lepuslapis
 *
 */

public class HelicesMaker {
	
		/**
		 * given an initial pair, looks around that pair to form a helix with all Deltas > 0
		 * @param leftIdx
		 * @param rightIdx
		 * @param pairedProbs
		 * @param unpairedProbs
		 * @param canPair
		 * @return said helix
		 */
		public Helix makeHelix(int leftIdx, int rightIdx, PointRes[][] diffs, boolean[][] canPair) {
			int length = diffs.length;	
			int helixLength = 1;
			while (diffs[leftIdx+helixLength][rightIdx-helixLength].signum()>0) {
				helixLength++;
			}
			while ((leftIdx>0)&&(rightIdx<length-1)&&(diffs[leftIdx-1][rightIdx+1].signum()>0)) {
				leftIdx--; rightIdx++; helixLength++;
			}
			Helix helix = new Helix(leftIdx,rightIdx,helixLength,diffs);
			return helix;
		}
		
		public Helix makeHelix(int leftIdx, int rightIdx, double[][] diffs, boolean[][] canPair) {
			int length = diffs.length;	
			int helixLength = 1;
			while (diffs[leftIdx+helixLength][rightIdx-helixLength]>0) {
				helixLength++;
			}
			while ((leftIdx>0)&&(rightIdx<length-1)&&(diffs[leftIdx-1][rightIdx+1]>0)) {
				leftIdx--; rightIdx++; helixLength++;
			}
			Helix helix = new Helix(leftIdx,rightIdx,helixLength,diffs);
			return helix;
		}
}

/*
 * PointRes[][] map = new PointRes[length][length];
			for (int j = 0; j<length; j++) {
				map[j][j] = unpairedProbs[j];
				for (int k = j+1; k<length; k++) {
					map[j][k]=PointRes.ZERO;
				}
			}
			boolean[][] mapsteps = new boolean[length][length];
			
			// start by doing MAP for subsequence leftIdx+1..rightIdx-1 
			for (int b = 1; b<rightIdx-leftIdx-1; b++) {
				for (int j = leftIdx + 1; j<rightIdx-b; j++) {
					map[j][j+b] = (unpairedProbs[j].divide(PointRes.valueOf(2))).add(map[j+1][j+b],MathContext.DECIMAL32);
					PointRes tmp = PointRes.ZERO;
					for (int k = j+1; k<j+b; k++) {
						tmp = pairedProbs[j][k].add(map[j+1][k-1]).add(map[k+1][j+b],MathContext.DECIMAL32);
						if (tmp.compareTo(map[j][j+b])>0) {
							map[j][j+b] = tmp; 
						}
					}
					tmp = pairedProbs[j][j+b].add(map[j+1][j+b-1]);
					if (tmp.compareTo(map[j][j+b])>0) {
						map[j][j+b] = tmp; mapsteps[j][j+b] = true;
					}
				}
			}
			int helixLength = 1;
			int j = 1;
			while (mapsteps[leftIdx+j][rightIdx-j]) {
				helixLength++;
			}
			
			//now do MAP for the other side of the potential helix; this is easier: we can work outwards progressively
			j = 1; int leftIdxtmp = leftIdx; int rightIdxtmp = rightIdx;
			PointRes leftup = PointRes.ZERO; PointRes rightup = PointRes.ZERO;
			PointRes pair = PointRes.ZERO;
			while ((leftIdx-j>=0)&&(rightIdx+j<length)) {
				leftup = (unpairedProbs[leftIdx-j].divide(PointRes.valueOf(2))).add(map[j+1][j+b],MathContext.DECIMAL32);
			}
			leftIdx = leftIdxtmp; rightIdx = rightIdxtmp;
*/