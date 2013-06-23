package uk.ac.ox.osscb.domain;


public class NucleotideProbsDouble extends SquareMatrixPlusVector<Double> {

	protected NucleotideProbsDouble(int dim) {
		super(dim, (double)0);// do not invoke method while directly invoking constructor
	}
	
	private char[]alphabet;
	
	public NucleotideProbsDouble(double[][] pairingProbs, double[] unpairingProbs, char[] alphabet){
		this(pairingProbs.length);
		
		if(null == unpairingProbs)
			throw new IllegalArgumentException("TODO-2");
		if(null == alphabet)
			throw new IllegalArgumentException("TODO-3");
		
		for(int rowIdx = 0; rowIdx < getDim(); rowIdx++){
			if(pairingProbs[rowIdx].length != getDim()){
				throw new IllegalArgumentException(String.format(
						"row %d (zbi), expected len to be %d, actual was %d. Must be same size",
						rowIdx, getDim(), pairingProbs[rowIdx].length));
			}
		}

		if(unpairingProbs.length != getDim()){
			throw new IllegalArgumentException(String.format(
					"Unpairing probs expected len to be %d, actual was %d. Must be same size",
					getDim(), unpairingProbs.length));
		}
		
		if(alphabet.length != getDim()){
			throw new IllegalArgumentException(String.format(
					"Alphabet expected len to be %d, actual was %d. Must be same size",
					alphabet, unpairingProbs.length));
		}
		
		for(int rowIdx = 0; rowIdx < getDim(); rowIdx++){
			for(int colIdx = 0; colIdx < getDim(); colIdx++){
				setPairingProbability(rowIdx, colIdx, pairingProbs[rowIdx][colIdx]);
			}
		}
		for(int rowIdx = 0; rowIdx < getDim(); rowIdx++){
			setUnpairingProbability(rowIdx, unpairingProbs[rowIdx]);
		}
		
		this.alphabet = new char[getDim()];
		System.arraycopy(alphabet, 0, this.alphabet, 0, getDim());
	}

	@Override
	protected Double getInitialNumber() {
		return (double)0;
	}

	public char[] getAlphabet() {
		// copy to avoid changing
		char[]res = new char[getDim()];
		System.arraycopy(this.alphabet, 0, res, 0, getDim());
		return res;
	}
	
	public String getAlphabetStr() {
		return new String(getAlphabet());
	}
}
