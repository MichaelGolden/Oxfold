package uk.ac.ox.osscb;

public class PossiblePairFinder {
	
	public boolean[][] canPair(int[] structure2, boolean [] delete) {
		if (null == structure2) {
			throw new IllegalArgumentException("Oh dear: the structure is null.");
		}
		
		int [] structure = CoFoldAnalogue.reinsertDeleted(structure2, delete, Constants.UnpairedBaseIdx);
		
		int length = structure.length;
		boolean[][] cp = new boolean[length][length];
		for (int j = 0; j<length-4; j++) {
			if (structure[j] == -1) {
				int k = j + 1;
				while (k < length) {
					if (structure[k] > k) {
						k = structure[k] + 1;
					} else if (structure[k] >= 0) {
						k = length;
					} else if (k > j + 3) {
						cp[j][k] = true; k++;
					} else {
						k++;
					}
				}
			} else if (structure[j]>j) {
				cp[j][structure[j]] = true;
			}
		}
		CoFoldAnalogue.deleteColumns(cp, delete);
		return cp;
	}
	/**
	 * finds possible pairs in structure
	 * @author lepuslapis
	 * @param structure is the structure matrix with s[j]=k if j and k are paired. 
	 * @return matrix cp with cp[j,k]=1 if j+4<=k and j and k can pair, and all other entries zero.
	 */
	public boolean[][] canPair(int[] structure) {
		if (null == structure) {
			throw new IllegalArgumentException("Oh dear: the structure is null.");
		}
		int length = structure.length;
		boolean[][] cp = new boolean[length][length];
		for (int j = 0; j<length-4; j++) {
			if (structure[j] == -1) {
				int k = j + 1;
				while (k < length) {
					if (structure[k] > k) {
						k = structure[k] + 1;
					} else if (structure[k] >= 0) {
						k = length;
					} else if (k > j + 3) {
						cp[j][k] = true; k++;
					} else {
						k++;
					}
				}
			} else if (structure[j]>j) {
				cp[j][structure[j]] = true;
			}
		}
		return cp;
	}
}
