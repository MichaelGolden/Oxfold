package uk.ac.ox.osscb.inoutside;

import java.math.BigDecimal;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;

/**
 * Inside outside implementation
 * 
 * @author Vladimir
 */
public class InsideOutsideCalculator implements IOsideCalculator {
	
	private Grammar grammar;

	public InsideOutsideCalculator(Grammar grammar) {
		super();

		if(null == grammar)
			throw new NullPointerException("grammar cannot be null");

		this.grammar = grammar;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#inside(uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, int[], boolean[][])
	 */
	public InsideOutsideProbabilities inside(NucleotideProbsPrecise pairingProbs, 
			int[][] distances, double weight, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					BigDecimal prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int h = j; h < j+b; h++) {
						BigDecimal prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						BigDecimal prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					BigDecimal probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double exp = Math.exp(-distances[j][j+b] / weight);
						BigDecimal probIncrement = BigDecimal.ZERO;
						probIncrement = rule2.getProbability().multiply(BigDecimal.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}			
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}

	public InsideOutsideProbabilities insideE(NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					BigDecimal prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int h = j; h < j+b; h++) {
						BigDecimal prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						BigDecimal prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					BigDecimal probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						BigDecimal probIncrement = BigDecimal.ZERO;
						probIncrement = rule2.getProbability().multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}			
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}
	
	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	public InsideOutsideProbabilities outside(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, int[][] distances, double weight, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, BigDecimal.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), j, k);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					BigDecimal probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = BigDecimal.ZERO;
					for(int k = 0; k<j;k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						BigDecimal probIncrement = BigDecimal.ZERO;
						BigDecimal exp = BigDecimal.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability().multiply(exp)
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
	}
	
	public InsideOutsideProbabilities outsideE(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, BigDecimal.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), j, k);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					BigDecimal probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = BigDecimal.ZERO;
					for(int k = 0; k<j;k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						BigDecimal probIncrement = BigDecimal.ZERO;
						probIncrement = rule2.getProbability().multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
	}

	@Override
	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					BigDecimal prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int h = j; h < j+b; h++) {
						BigDecimal prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						BigDecimal prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					BigDecimal probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double exp = Math.exp(-distances[j][j+b]/ weight);
						BigDecimal probIncrement = BigDecimal.ZERO;
						probIncrement = rule2.getProbability().multiply(BigDecimal.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}			
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}

	@Override
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
		
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, BigDecimal.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					BigDecimal tmp = BigDecimal.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), j, k);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					BigDecimal probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = BigDecimal.ZERO;
					for(int k = 0; k<j;k++) {
						BigDecimal oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						BigDecimal iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						BigDecimal probIncrement = BigDecimal.ZERO;
						BigDecimal exp = BigDecimal.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability().multiply(exp)
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
		
	}
}
