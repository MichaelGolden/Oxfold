package uk.ac.ox.osscb.domain;

import java.math.BigDecimal;

import uk.ac.ox.osscb.Util;

/**
 * Diagonal matrix + unpairing probabilities
 * 
 * @author Vladimir
 */
@Deprecated
public class NucleotideProbsPreciseOld {
	
	private int dim;
	
	private BigDecimal[][] probs;
	private BigDecimal[] unpairedProb;
	
	public NucleotideProbsPreciseOld(int dim) {
		super();
		if(dim < 1){
			throw new RuntimeException("");
		}
		
		this.dim = dim;
		this.probs = Util.makeSquareZerosMatrix(dim);
		
		this.unpairedProb = new BigDecimal[dim];
		for(int dimIdx = 0; dimIdx < dim; dimIdx++){
			this.unpairedProb[dimIdx] = BigDecimal.valueOf(0);
		}
	}

	public BigDecimal setPairingProbability(int i, int j, double prob){
		return setPairingProbability(i, j, BigDecimal.valueOf(prob));
	}

	public BigDecimal setPairingProbability(int i, int j, BigDecimal prob){
		checkDimensions(i, j);
		
		// TODO: make the mtx be triangular for memory efficiency 
		// (do we care since we've got several hundreds of nucleotides typically?)
		// set them both, will take only one later on.
		this.probs[i][j] = prob;
		this.probs[j][i] = prob;
		return prob;
	}

	public BigDecimal getPairingProbability(int i, int j){
		checkDimensions(i, j);
		return this.probs[i][j];
	}
	
	public BigDecimal getUnpairingProbability(int i){
		checkDimensions(i);
		return this.unpairedProb[i];
	}
	
	public BigDecimal setUnpairingProbability(int i, double prob){
		return setUnpairingProbability(i, BigDecimal.valueOf(prob));
	}
	
	public BigDecimal setUnpairingProbability(int i, BigDecimal prob){
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
