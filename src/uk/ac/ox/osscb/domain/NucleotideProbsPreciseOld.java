package uk.ac.ox.osscb.domain;



import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.Util;

/**
 * Diagonal matrix + unpairing probabilities
 * 
 * @author Vladimir
 */
@Deprecated
public class NucleotideProbsPreciseOld {
	
	private int dim;
	
	private PointRes[][] probs;
	private PointRes[] unpairedProb;
	
	public NucleotideProbsPreciseOld(int dim) {
		super();
		if(dim < 1){
			throw new RuntimeException("");
		}
		
		this.dim = dim;
		this.probs = Util.makeSquareZerosMatrix(dim);
		
		this.unpairedProb = new PointRes[dim];
		for(int dimIdx = 0; dimIdx < dim; dimIdx++){
			this.unpairedProb[dimIdx] = PointRes.valueOf(0);
		}
	}

	public PointRes setPairingProbability(int i, int j, double prob){
		return setPairingProbability(i, j, PointRes.valueOf(prob));
	}

	public PointRes setPairingProbability(int i, int j, PointRes prob){
		checkDimensions(i, j);
		
		// TODO: make the mtx be triangular for memory efficiency 
		// (do we care since we've got several hundreds of nucleotides typically?)
		// set them both, will take only one later on.
		this.probs[i][j] = prob;
		this.probs[j][i] = prob;
		return prob;
	}

	public PointRes getPairingProbability(int i, int j){
		checkDimensions(i, j);
		return this.probs[i][j];
	}
	
	public PointRes getUnpairingProbability(int i){
		checkDimensions(i);
		return this.unpairedProb[i];
	}
	
	public PointRes setUnpairingProbability(int i, double prob){
		return setUnpairingProbability(i, PointRes.valueOf(prob));
	}
	
	public PointRes setUnpairingProbability(int i, PointRes prob){
		checkDimensions(i);
		this.unpairedProb[i] = prob;
		return prob;
	}
	
	private void checkDimensions(int i) {
		if(i >= this.dim){
			throw new RuntimeException("TODO-3");
		}
	}

	private void checkDimensions(int i, int j) {
		if(i >= this.dim){
			throw new RuntimeException("TODO-1");
		}
		if(j >= this.dim){
			throw new RuntimeException("TODO-2");
		}
	}
}
