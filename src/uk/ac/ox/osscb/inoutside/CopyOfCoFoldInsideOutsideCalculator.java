package uk.ac.ox.osscb.inoutside;



import java.io.File;
import java.io.IOException;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;

/**
 * Inside outside implementation
 * 
 * @author Vladimir
 */
public class CopyOfCoFoldInsideOutsideCalculator {
	
	private Grammar grammar;

	public CopyOfCoFoldInsideOutsideCalculator(Grammar grammar) {
		super();

		if(null == grammar)
			throw new NullPointerException("grammar cannot be null");

		this.grammar = grammar;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#inside(uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, int[], boolean[][])
	 */
	public InsideOutsideProbabilities inside(NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						try {
							IO.writeLine(new File("inside.txt"), "R3\t"+j+"\t"+h, true, true);
							IO.writeLine(new File("inside.txt"), "R3\t"+(h+1)+"\t"+(j+b), true, true);
						} catch (IOException e) {
							e.printStackTrace();
						}
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						try {
							IO.writeLine(new File("inside.txt"), "R2\t"+j+"\t"+(j+b), true, true);
							IO.writeLine(new File("inside.txt"), "R2\t"+(j+1)+"\t"+(j+b-1), true, true);
						} catch (IOException e) {
							e.printStackTrace();
						}
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
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
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
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
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						double distance = Math.abs(b);
						PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
						//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
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
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
	}

	//@Override
	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double alpha, double tau, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
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

	//@Override
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
		
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						double distance = Math.abs(b);
						PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
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
