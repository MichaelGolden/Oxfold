package uk.ac.ox.osscb.inoutside;


import java.math.RoundingMode;
import java.util.HashMap;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;

public class PosteriorProbabilitiesCalculator {
		
		private Grammar grammar;
		
		public PosteriorProbabilitiesCalculator(Grammar grammar) {
			this.grammar = grammar;
		}
		
		public PosteriorProbabilities calculate(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				int[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							pairedProbs[j][k] = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);

							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}


		public ShortDoublePosteriorProbabilities calculateShort(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			double[][] pairedProbs = new double[length][length];
			double[] unpairedProbs = new double[length];
			
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					PointRes prob = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
					unpairedProbs[j] = prob.doubleValue();
				}
			}
			//paired probabilities
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							PointRes prob = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
							pairedProbs[j][k] = prob.doubleValue();
						}
					}
				}
			}
			ShortDoublePosteriorProbabilities posteriorProbabilities = new ShortDoublePosteriorProbabilities(unpairedProbs, pairedProbs);
			return posteriorProbabilities;
		}
		
		public PosteriorProbabilities calculateE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							pairedProbs[j][k] = tmp.multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
													
							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}
		
		//takes double[][] for distance param (not int[][])
		public PosteriorProbabilities calculate(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							pairedProbs[j][k] = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);

							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}
		
		public ShortDoublePosteriorProbabilities calculateShortE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			double[][] pairedProbs = new double[length][length];
			double[] unpairedProbs = new double[length];
			
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					PointRes prob = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
					unpairedProbs[j] = prob.doubleValue();
				}
			}
			//paired probabilities
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes prob = tmp.multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
							pairedProbs[j][k] = prob.doubleValue();
						}
					}
				}
			}
			ShortDoublePosteriorProbabilities posteriorProbabilities = new ShortDoublePosteriorProbabilities(unpairedProbs, pairedProbs);
			return posteriorProbabilities;
		}
		
		public PointRes calculateCompE(boolean[][] incomp, Helix helix, InsideOutsideProbabilities insideProbs, InsideOutsideProbabilities outsideProbs,
				NucleotideProbsPrecise nucleotideProbs, PointRes[][] pairedProbs) {
			PointRes rprob = PointRes.ZERO;
			for (int j=0; j<incomp.length; j++) {
				for (int k = j+1; k<incomp.length; k++) {
					if (incomp[j][k]) {
						rprob = rprob.add(pairedProbs[j][k]);
					}
				}
			}
			PointRes hprob = findHelixProbabilityE(insideProbs, outsideProbs, nucleotideProbs, helix);
			if (rprob.signum()>0) {
				rprob = hprob.divide(rprob.add(hprob), RoundingMode.HALF_UP);
			}
			return rprob;
		}
		
		/**
		 * calculate accuracy scores
		 */
		public PointRes[][] getDiffs(PointRes[][] pairedProbs, PointRes[] unpairedProbs, boolean[][] canPair) {
			int length = unpairedProbs.length;
			PointRes[][] diffs = new PointRes[length][length];
			for (int j = 0; j<length; j++) {
				for (int k = 0; k<length; k++) {
					if (canPair[j][k]) {
						diffs[j][k] = pairedProbs[j][k].subtract((unpairedProbs[j].add(unpairedProbs[k])).divide(PointRes.valueOf(2)));
					} else {
						diffs[j][k] = PointRes.valueOf(-1);
					}	
				}
			}
			return diffs;
		}
		
		
		/**
		 * calculate probability of helix forming
		 */
		
		private PointRes findHelixProbabilityE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, Helix helix) {
			int left = helix.getLeftIdx(); int right = helix.getRightIdx();
			int length = helix.getHelixLength();
			PointRes nprob = PointRes.ONE;
			for (int j = 0; j<length; j++) {
				nprob = nprob.multiply(nucleotideProbs.getPairingProbability(left+j, right-j));
			}
			nprob = nprob.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, insideProbs.getDimension()-1), RoundingMode.HALF_UP);
			HashMap<Character,PointRes> Probs = new HashMap<Character,PointRes>();
			//initialise Hashmap with inside probabilities
			for (ProductionRule pr : this.grammar.getRules(RuleType.RULE2)) {
				Probs.put(pr.getLeft(), insideProbs.getProb(pr.getRight()[1], left+length, right-length).multiply(pr.getProbability())); 
			}
			for (int j = length-2; j>=0; j--) {
				HashMap<Character,PointRes> tmpProbs = new HashMap<Character,PointRes>();
				for (char nt: this.grammar.getNonterminals()) {
					tmpProbs.put(nt, PointRes.ZERO);
				}
				for (ProductionRule pr: this.grammar.getRules(RuleType.RULE2)) {
					PointRes newtmp = pr.getProbability().multiply(Probs.get(pr.getRight()[1]));
					tmpProbs.put(pr.getLeft(), newtmp);
				}
				Probs = tmpProbs;
			}
			PointRes prob = PointRes.ZERO;
			for (char nt : this.grammar.getNonterminals()) {
				prob = prob.add(outsideProbs.getProb(nt, left, right).multiply(Probs.get(nt)));
			}
			return prob;
		}
}
