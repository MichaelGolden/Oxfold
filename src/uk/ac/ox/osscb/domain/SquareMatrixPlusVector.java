package uk.ac.ox.osscb.domain;

import java.lang.reflect.Array;


/**
 * Generic class incapsulating a triangle matrix. In order not to bother
 * with indices convertion we keep a[i,j] == a[j,i]
 * 
 * @author Vladimir
 *
 * @param <NT> type of number we are using. E.g. are {@link PointRes},
 * {@link Double}
 */
public abstract class SquareMatrixPlusVector<NT extends Number> {

	private int dim;
	
	private NT[][] probs;
	private NT[] unpairedProb;
	
	/**
	 * 
	 * @param dim
	 * @param exampleValue used to derive runtime {@link Class<NT>} from.
	 */
	@SuppressWarnings("unchecked")
	protected SquareMatrixPlusVector(int dim, NT exampleValue) {
		super();
		if(dim < 1){
			throw new RuntimeException("");
		}
		
		this.dim = dim;
		
		this.unpairedProb = (NT[])Array.newInstance(exampleValue.getClass(), dim);
		for(int dimIdx = 0; dimIdx < dim; dimIdx++){
			this.unpairedProb[dimIdx] = getInitialNumber();
		}
		
		this.probs = (NT[][])createNumberArray(dim);
		
		for(int i = 0; i < dim; i++){
			for(int j = 0; j < dim; j++){
				// Object numberArray = createNumberArray(dim);
				
				this.probs[i] = (NT[])Array.newInstance(exampleValue.getClass(), dim);
				this.probs[i][j] = getInitialNumber();
			}
		}
	}

	private Object createNumberArray(int dim) {
		return Array.newInstance(this.unpairedProb.getClass(), dim);
	}
	
	protected abstract NT getInitialNumber();
	
	public NT setUnpairingProbability(int i, NT prob){
		checkDimensions(i);
		this.unpairedProb[i] = prob;
		return prob;
	}
	
	public NT setPairingProbability(int i, int j, NT prob){
		checkDimensions(i, j);
		
		this.probs[i][j] = prob;
		return prob;
	}

	public NT[][] getMtx(){
		return this.probs;
	}

	public NT[] getVector(){
		return this.unpairedProb;
	}

	public NT getPairingProbability(int i, int j){
		checkDimensions(i, j);
		return this.probs[i][j];
	}
	
	public NT getUnpairingProbability(int i){
		checkDimensions(i);
		return this.unpairedProb[i];
	}
	
	/* public NT getUnpairedColumnProduct(int j){
	checkDimensions(j, -1);

	NT product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
}

public NT getPairedColumnProduct(int j, int k){
	checkDimensions(j, k);
	NT product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
}*/	

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

	public int getDim() {
		return dim;
	}
}
