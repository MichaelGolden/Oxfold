package uk.ac.ox.osscb;

/**
 * Class for the data object in nodes in an evolutionary tree.
 * @author lmath
 *
 */
public class EvolutionaryNodeData {
	private double branchLength;
	private int idx;
	
	public EvolutionaryNodeData(double branchLength)
	{
		setBranchLength(branchLength);
		idx = -1;
	}
	
	public EvolutionaryNodeData()
	{
		branchLength = -1;
		idx = -1;
	}
	

	public int getIdx() {
		
		return this.idx;
	}

	public void setIdx(int idx) {
	//	if(-1 != this.idx)
	//		throw new IllegalStateException(String.format("Tried to reset an index " +
	//				"in Evolutionary Tree Traversal to %d. " +
	//				"Index is already set to %d", idx, this.idx));
		this.idx = idx;
	}

	public double getBranchLength(){
		
		return this.branchLength;
	}
	
	public void setBranchLength(double branchLength){
		
		this.branchLength = branchLength;
	}
	
}
