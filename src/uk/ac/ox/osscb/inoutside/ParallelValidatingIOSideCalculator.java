package uk.ac.ox.osscb.inoutside;


//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;

public class ParallelValidatingIOSideCalculator implements IOsideCalculator{
	
	//private final Logger log = LoggerFactory.getLogger(ValidatingIOSideCalculator.class);
	
	private ParallelInsideOutsideCalculator ioCalcInternal;
	private IOsideCalculator ioCalcInternal2;
	private Grammar grammar;

	public ParallelValidatingIOSideCalculator(ParallelInsideOutsideCalculator iOsideCalculator, IOsideCalculator iOsideCalculator2, Grammar grammar) {
		super();
		this.ioCalcInternal = iOsideCalculator;
		this.ioCalcInternal2 = iOsideCalculator2;
		this.grammar = grammar;
	}

	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		
		InsideOutsideProbabilities inside = this.ioCalcInternal.insideParallelOxfold(pairingProbs, distances, weight,
						structure, canPair);
		
		return validateInside(inside);
	}
	
	public InsideOutsideProbabilities insideE(
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		
		InsideOutsideProbabilities inside = this.ioCalcInternal.insideParallel(pairingProbs, structure, canPair, null);
		//InsideOutsideProbabilities inside = this.ioCalcInternal2.insideE(pairingProbs, structure, canPair);
		
		return validateInside(inside);
	}

	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outsideParallelOxfold(insideProbs, pairingProbs, distances,
						weight, structure, canPair);
		return validateOutside(outside);
	}

	public InsideOutsideProbabilities outsideE(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outsideParallel(insideProbs, pairingProbs, structure, canPair, null);
		//InsideOutsideProbabilities outside = this.ioCalcInternal2.outsideE(insideProbs, pairingProbs, structure, canPair);
		return validateOutside(outside);
	}
	
	private InsideOutsideProbabilities validateInside(InsideOutsideProbabilities inside) {
		
		PointRes divisor = inside.getProb(this.grammar.getNonterminals()[0], 0, inside.getDimension()-1);
		
		if(0 == divisor.compareTo(PointRes.ZERO) ){
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
		InsideOutsideProbabilities inside = this.ioCalcInternal.insideParallelOxfold(pairingProbs, distances, weight,
				structure, canPair);

		return validateInside(inside);

	}

	
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		InsideOutsideProbabilities outside = this.ioCalcInternal.outsideParallelOxfold(insideProbs, pairingProbs, distances,
				weight, structure, canPair);
		
		return validateOutside(outside);
	}
}
