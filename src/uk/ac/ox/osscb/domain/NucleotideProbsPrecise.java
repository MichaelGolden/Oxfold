package uk.ac.ox.osscb.domain;

import java.math.BigDecimal;
import java.math.MathContext;

import uk.ac.ox.osscb.util.ProbabilityValueValidator;

public class NucleotideProbsPrecise extends SquareMatrixPlusVector<BigDecimal> {
	
	private MathContext mathCtx = null;//new MathContext(5);


	public NucleotideProbsPrecise(int dim) {
		super(dim, BigDecimal.ZERO);
	}

	@Override
	protected BigDecimal getInitialNumber() {
		return BigDecimal.ZERO;
	}

	@Override
	public BigDecimal setUnpairingProbability(int i, BigDecimal prob) {
		ProbabilityValueValidator.validateP(prob);
		if(null != this.mathCtx){
			prob = prob.round(this.mathCtx);
		}
		return super.setUnpairingProbability(i, prob);
	}

	@Override
	public BigDecimal setPairingProbability(int i, int j, BigDecimal prob) {
		ProbabilityValueValidator.validateP(prob);
		if(null != this.mathCtx){
			prob = prob.round(this.mathCtx);
		}
		return super.setPairingProbability(i, j, prob);
	}
	
	
	
	/* 
	public BigDecimal getUnpairedColumnProduct(int j){
	checkDimensions(j, -1);

	BigDecimal product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
}

public BigDecimal getPairedColumnProduct(int j, int k){
	checkDimensions(j, k);
	BigDecimal product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
	}
*/
	
}
