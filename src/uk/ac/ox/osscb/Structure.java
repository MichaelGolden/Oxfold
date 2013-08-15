package uk.ac.ox.osscb;

import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.PPOutput;

/*
 * Angela messing around...
 */
public class Structure {

	private int [] pairings; 
	private boolean [] keepPairs; 
	
	public Structure(int[] struct){
		pairings = struct; 
		keepPairs = new boolean[pairings.length];
	}
	
	public Structure(int[] struct, boolean[] b){
		pairings = struct; 
		keepPairs = b;
	}
	
	public int[] makeNewStructure(PPOutput ppProbs, boolean keepFlag){
		for (int j=0;j<ppProbs.gethelixLength(); j++) {
			pairings[ppProbs.getLeftIdx()+j] = ppProbs.getRightIdx()-j;
			pairings[ppProbs.getRightIdx()-j] = ppProbs.getLeftIdx()+j;
			keepPairs[ppProbs.getLeftIdx()+j] = keepFlag; 
			keepPairs[ppProbs.getRightIdx()-j] = keepFlag;
		}
		return pairings;
	}
	
	public int[] makeNewStructure(Helix helix, boolean keepFlag) {
		for (int j=0;j<helix.getHelixLength(); j++) {
			pairings[helix.getLeftIdx()+j] = helix.getRightIdx()-j;
			pairings[helix.getRightIdx()-j] = helix.getLeftIdx()+j;
			keepPairs[helix.getLeftIdx()+j] = keepFlag; 
			keepPairs[helix.getRightIdx()-j] = keepFlag;
		}
		return pairings;
	}
	
	public void setPairings(int[] tmp){
		pairings = tmp; 
	}
	
	public void updateKeepPairs(int i, boolean f){
		keepPairs[i] = f;
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
	
	public int[] getPairings(){
		return pairings; 
	}
	
	public boolean[] getKeepPairs(){
		return keepPairs; 
	}
	
	public int[] removeWeakPairs(){
		for(int i = 0; i < pairings.length; i++){
			if(keepPairs[i]==false){
				pairings[i] = -1; 
			}
		}
		return pairings; 
	}
}
