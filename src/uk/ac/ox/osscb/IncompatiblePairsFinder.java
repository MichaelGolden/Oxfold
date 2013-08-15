package uk.ac.ox.osscb;


import java.math.RoundingMode;
import java.util.LinkedList;

import uk.ac.ox.osscb.inoutside.Helix;

/**
 * 
 * @author Lepuslapis
 *
 */

public class IncompatiblePairsFinder {
	
		/**
		 * determines lower bound on compatibility of pair left right
		 * @param canPair
		 * @param left
		 * @param right
		 * @param pairedProbs
		 * @return lower bound on compatibility
		 */
		public PointRes calculateComp(boolean[][] incomp, int left, int right, PointRes[][] pairedProbs) {
			PointRes rprob = PointRes.ZERO;
			for (int j=0; j<incomp.length; j++) {
				for (int k = j+1; k<incomp.length; k++) {
					if (incomp[j][k]) {
						rprob = rprob.add(pairedProbs[j][k]);
					}
				}
			}
			if (rprob.signum()>0) {
				rprob = pairedProbs[left][right].divide(rprob, RoundingMode.HALF_UP);
			}
			return rprob;
		}
	
		/**
		 * determines possible pairs that are incompatible with the pair left right
		 * @param canPair
		 * @param left
		 * @param right
		 * @return boolean array of incompatible pairs
		 */
		
		public boolean[][] find(boolean[][] canPair, int left, int right) {
			boolean[][] incomp = new boolean[canPair.length][canPair.length];
			for (int j = 0; j<left; j++) {
				for (int k = left; k<=right; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
			}
			for (int k = left+1; k<canPair.length; k++) {
				if (canPair[left][k]) {
					incomp[left][k] = true;
				}
			}
			for (int j = left+1; j<=right; j++) {
				for (int k = right; k<canPair.length; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
			}
			return incomp;
		}
		
		public boolean[][] findAll(boolean[][] canPair, Helix helix) {
			int left = helix.getLeftIdx(); int right = helix.getRightIdx();
			int length = helix.getHelixLength();
			boolean[][] incomp = new boolean[canPair.length][canPair.length];
			for (int j = 0; j<left; j++) {
				for (int k = left; k<=right; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
			}
			for (int j = left; j<left+length; j++) {
				for (int k = j+1; k<right-j+left; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
				for (int k = right-j+left+2; k<canPair.length; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
			}
			for (int j = left+length; j<=right; j++) {
				for (int k = Math.max(right-length+1, j+1); k<canPair.length; k++) {
					if (canPair[j][k]) {
						incomp[j][k] = true;
					}
				}
			}
			return incomp;
		}
		
		/*public LinkedList<Helix> getIncompatibleHelices(boolean[][] canPair, boolean[][] incomp, PointRes[][] diffs, 
				int leftIdx, int rightIdx) {
			LinkedList<Helix> Helices = new LinkedList<Helix>();
			boolean[][] tmpincomp = new boolean[incomp.length][incomp.length];
			System.arraycopy(incomp, 0, tmpincomp, 0, incomp.length);
			int length = canPair.length;
			for (int j = 0; j<length; j++) {
				for (int k = j+1; k<length; k++) {
					if (((j != leftIdx)||(k != rightIdx))&&(tmpincomp[j][k])&&(diffs[j][k].signum()>0)) {
						Helix newHelix = new HelicesMaker().makeHelix(j, k, diffs, canPair);
						if ((newHelix.getHelixLength()>=1)) {
							for (int l = 0; l<newHelix.getHelixLength(); l++) {
								tmpincomp[newHelix.getLeftIdx()+l][newHelix.getRightIdx()-l] = false;
							}
							Helices.add(newHelix);
						}
					}
				}
			}
			return Helices;
		}*/
		public LinkedList<Helix> getIncompatibleHelices(boolean[][] canPair, boolean[][] incomp, PointRes [][] pairedProbs, PointRes[][] diffs, 
				int leftIdx, int rightIdx) {
			LinkedList<Helix> Helices = new LinkedList<Helix>();
			boolean[][] tmpincomp = new boolean[incomp.length][incomp.length];
			System.arraycopy(incomp, 0, tmpincomp, 0, incomp.length);
			int length = canPair.length;
			for (int j = 0; j<length; j++) {
				for (int k = j+1; k<length; k++) {
					if (((j != leftIdx)||(k != rightIdx))&&(tmpincomp[j][k])&&(diffs[j][k].signum()>0)) {
						//Helix newHelix = new HelicesMaker().makeHelix(j, k, pairedProbs, diffs, canPair);
						Helix newHelix = new HelicesMaker().makeHelix(j, k, diffs, canPair);
						if ((newHelix.getHelixLength()>=1)) {
							for (int l = 0; l<newHelix.getHelixLength(); l++) {
								tmpincomp[newHelix.getLeftIdx()+l][newHelix.getRightIdx()-l] = false;
							}
							Helices.add(newHelix);
						}
					}
				}
			}
			return Helices;
		}
		
}
