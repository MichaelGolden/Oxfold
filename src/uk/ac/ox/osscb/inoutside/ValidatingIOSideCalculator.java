package uk.ac.ox.osscb.inoutside;

import java.math.BigDecimal;

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;

public class ValidatingIOSideCalculator implements IOsideCalculator{
	
	//private final Logger log = LoggerFactory.getLogger(ValidatingIOSideCalculator.class);
	
	private IOsideCalculator ioCalcInternal;
	private Grammar grammar;

	public ValidatingIOSideCalculator(IOsideCalculator iOsideCalculator, Grammar grammar) {
		super();
		this.ioCalcInternal = iOsideCalculator;
		this.grammar = grammar;
	}

	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		
		InsideOutsideProbabilities inside = this.ioCalcInternal.inside(pairingProbs, distances, weight,
						structure, canPair);
		
		return validateInside(inside);
	}
	
	public InsideOutsideProbabilities insideE(
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		
		InsideOutsideProbabilities inside = this.ioCalcInternal.insideE(pairingProbs, structure, canPair);
		
		return validateInside(inside);
	}

	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outside(insideProbs, pairingProbs, distances,
						weight, structure, canPair);
		return validateOutside(outside);
	}

	public InsideOutsideProbabilities outsideE(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outsideE(insideProbs, pairingProbs, structure, canPair);
		return validateOutside(outside);
	}
	
	private InsideOutsideProbabilities validateInside(InsideOutsideProbabilities inside) {
		
		BigDecimal divisor = inside.getProb(this.grammar.getNonterminals()[0], 0, inside.getDimension()-1);
		
		if(0 == divisor.compareTo(BigDecimal.ZERO) ){
			throw new UnsupportedOperationException(String.format("About to return zero divisor: %.30g", divisor));			
		}

		return inside;
	}

	private InsideOutsideProbabilities validateOutside(InsideOutsideProbabilities outside) {
		return outside;
	}


	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities inside = this.ioCalcInternal.inside(pairingProbs, distances, weight,
				structure, canPair);

		return validateInside(inside);

	}

	
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outside(insideProbs, pairingProbs, distances,
				weight, structure, canPair);
		
		return validateOutside(outside);
	}
}
