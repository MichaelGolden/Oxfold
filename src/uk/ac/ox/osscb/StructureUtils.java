package uk.ac.ox.osscb;

import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.PPOutput;

public class StructureUtils {
		
		/**
		 * update structure by including given helix
		 * @param oldStructure
		 * @param ppProbs or helix
		 * @return updated structure
		 */
		public int[] makeNewStructure(int[] oldStructure, PPOutput ppProbs) {
			int[] newStructure = new int[oldStructure.length];
			System.arraycopy(oldStructure, 0, newStructure, 0, oldStructure.length);
			for (int j=0;j<ppProbs.gethelixLength(); j++) {
				newStructure[ppProbs.getLeftIdx()+j] = ppProbs.getRightIdx()-j;
				newStructure[ppProbs.getRightIdx()-j] = ppProbs.getLeftIdx()+j;
			}
			return newStructure;
		}
		
		public int[] makeNewStructure(int[] oldStructure, Helix helix) {
			int[] newStructure = new int[oldStructure.length];
			System.arraycopy(oldStructure, 0, newStructure, 0, oldStructure.length);
			for (int j=0;j<helix.getHelixLength(); j++) {
				newStructure[helix.getLeftIdx()+j] = helix.getRightIdx()-j;
				newStructure[helix.getRightIdx()-j] = helix.getLeftIdx()+j;
			}
			return newStructure;
		}
		
		/**
		 * inverse algorithm to dumpstructure
		 */
		
		public int[] getStructureFromString(String strucString) {
			int[] structure = new int[strucString.length()];
			int last = 0; int lasttmp = 0;
			for (int j = 0; j<strucString.length(); j++) {
				char s = strucString.charAt(j);
				if (s == '.') {
					structure[j] = -1;
				} else if (s == '(') {
					structure[j] = last;
					last = j;
				} else if (s == ')') {
					lasttmp = structure[last];
					structure[last] = j;
					structure[j] = last;
					last = lasttmp;
				} else {
					throw new IllegalArgumentException("Illegal character in structure string.");
				}
			}
			return structure;
		}
		
		/**
		 * get substructure of given structure
		 * @param structure
		 * @param lPos
		 * @param rPos
		 * @return substructure array
		 */
		public int[] getSubStructure(int[] structure, int lPos, int rPos) {
			if (rPos > structure.length -1) {
				rPos = structure.length -1;
			}
			if (lPos>=rPos) {
				throw new IllegalArgumentException("Illegal Arguments for Substructure.");
			}
			int[] subStruc = new int[rPos-lPos+1];
			for (int j = lPos; j<=rPos; j++) {
				if ((structure[j]<=rPos)&&(structure[j]>=lPos)) {
					subStruc[j-lPos] = structure[j]-lPos;
				} else {
					subStruc[j-lPos] = -1;
				}
			}
			return subStruc;
		}
		
}
