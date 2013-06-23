package uk.ac.ox.osscb.inoutside;

import java.math.BigDecimal;
import java.util.LinkedList;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.DistancesCalculator;
import uk.ac.ox.osscb.HelicesMaker;
import uk.ac.ox.osscb.IncompatiblePairsFinder;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.LoggingOutputGenerator;
import uk.ac.ox.osscb.PossiblePairFinder;
import uk.ac.ox.osscb.StructureUtils;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;

public class DynamicPPProbabilitiesCalculator {
	
	//private final Logger log = LoggerFactory.getLogger(PPProbabilitiesCalculator.class);

	private Grammar grammar;
	
	public DynamicPPProbabilitiesCalculator(Grammar grammar) {
		this.grammar = grammar;
	}
	
	/**
	 * 
	 * @param insideProbs - 'i', output of {@link InsideOutsideCalculator#inside(NucleotideProbsPrecise, int[][], double, int[], boolean[][])}
	 * @param outsideProbs - 'o', output of {@link InsideOutsideCalculator#outside(InsideOutsideProbabilities, NucleotideProbsPrecise, int[][], double, boolean[][])}
	 * @param nucleotideProbs - 'n', the same as input to  {@link InsideOutsideCalculator#inside(NucleotideProbsPrecise, int[][], double)}
	 * @param distances - 'd', the same as input to  {@link InsideOutsideCalculator#inside(NucleotideProbsPrecise, int[][], double)}
	 * @param weight - 
	 * @param structure - the pairing structure to be calculated
	 * @param canPair - output of canpair method
	 */

	
	public PPOutput calculate(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs,
			int[][] distances, double weight, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculate(insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		BigDecimal[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		BigDecimal[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutput output = null;
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutput(-1,-1,0,BigDecimal.ZERO,BigDecimal.ZERO);
		} else {
			BigDecimal[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			BigDecimal diff = diffs[leftIdx][rightIdx];
			Helix mainHelix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			BigDecimal rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			Helix chosenHelix = mainHelix;
			int mainRight = mainHelix.getRightIdx();
			BigDecimal mainScore = mainHelix.getScore();
			// test whether we should be more careful about our choice of pairing...
			if (rprob.compareTo(Constants.DynamicCutOff)<0) {
				incomp = new IncompatiblePairsFinder().findAll(canPair, mainHelix);
				System.out.println("Possible incompatibility (" + rprob.doubleValue() + ") detected... Proposed new structure is (score " + mainScore.doubleValue() +"):");
				int[] tmpstructure = new StructureUtils().makeNewStructure(structure, chosenHelix);
				System.out.println(new LoggingOutputGenerator().dumpStructure(tmpstructure));
				LinkedList<Helix> incompHelices = new IncompatiblePairsFinder().getIncompatibleHelices(canPair,incomp,diffs,leftIdx,rightIdx);
				LinkedList<Helix> interestingHelices = new LinkedList<Helix>();
				System.out.println("Have determined " + incompHelices.size() + " incompatible helices.");
				for (Helix competingHelix: incompHelices) {
					System.out.println("Testing alternative structure: (score " + competingHelix.getScore().doubleValue() +"):");
					tmpstructure = new StructureUtils().makeNewStructure(structure, competingHelix);
					System.out.println(new LoggingOutputGenerator().dumpStructure(tmpstructure));
					int compRight = competingHelix.getRightIdx();
					if (competingHelix.getScore().compareTo(mainScore)>0) { //TODO: replace this with some threshold
						if (mainRight < compRight) {
							System.out.println("Competing structure scores more highly and competing helix further downstream. Breaking proposed helix.");
							// "break" main helix and form competing Helix
							interestingHelices.add(competingHelix);
						} else {
							System.out.println("Competing structure upstream from proposed structure, but scores more highly. Testing for transient structures.");
							// see whether we can prevent the competing Helix from forming
							int compLeft = competingHelix.getLeftIdx();
							int[] subStruc = new StructureUtils().getSubStructure(structure,0,compRight+Constants.TranscriptionOffset); 
							// now get Helices competing with Substructure
							SubOutput subOutput = getSubCompHelices(subStruc, nucleotideProbs, weight, compLeft, compRight);
							Helix subMainHelix = subOutput.getSubMainHelix();
							LinkedList<Helix> subCompHelices = subOutput.getSubCompHelices();
							BigDecimal subMainScore = subMainHelix.getScore();
							boolean found = false;
							for (Helix subCompHelix: subCompHelices) {
								if (subCompHelix.getScore().compareTo(subMainScore)>0) {
									found = true; break;
								}
							}
							if (!found) {
								interestingHelices.add(competingHelix);
								System.out.println("No transient structure found. Breaking main helix.");
							} else {
								System.out.println("Found transient structure. Breaking competing helix.");
							}
						}
					} else if (mainRight < compRight) {
						System.out.println("Proposed helix scores more highly than competing helix downstream. Testing for transient structures.");
						// see whether we can prevent the main Helix from forming
						int mainLeft = mainHelix.getLeftIdx();
						int[] subStruc = new StructureUtils().getSubStructure(structure,0,mainRight+Constants.TranscriptionOffset); 
						// now get matrix of probabilities for substructure
						SubOutput subOutput = getSubCompHelices(subStruc, nucleotideProbs, weight, mainLeft, mainRight);
						Helix subMainHelix = subOutput.getSubMainHelix();
						LinkedList<Helix> subCompHelices = subOutput.getSubCompHelices();
						BigDecimal subMainScore = subMainHelix.getScore();
						boolean found = false;
						for (Helix subCompHelix: subCompHelices) {
							if (subCompHelix.getScore().compareTo(subMainScore)>0) {
								found = true; break;
							}
						}
						if (found) {
							interestingHelices.add(competingHelix);
							System.out.println("Transient structure found. Breaking proposed helix.");
						} else {
							System.out.println("No transient structure found. Sticking with proposed helix.");
						}
						// otherwise "break" competing helix, i.e. do not include competing helix in list of interesting helices
					}
				}	
				// if no interesting structures: stick with main structure; otherwise compare interesting competing structures and decide.
				if (!interestingHelices.isEmpty()) {
					BigDecimal maxScore = BigDecimal.ZERO;
					chosenHelix = mainHelix;
					for (Helix helix: interestingHelices) {
						if (maxScore.compareTo(helix.getScore())<0) {
							chosenHelix = helix;
							maxScore = helix.getScore();
						}
					}	
				} else {
					chosenHelix = mainHelix;
				}
			} else {
				chosenHelix = mainHelix;
			}
			output = new PPOutput(chosenHelix.getLeftIdx(), chosenHelix.getRightIdx(), chosenHelix.getHelixLength(), diff, rprob);			
		}
		return output;
	}			
	
	private class SubOutput {
		
		private final Helix subMainHelix;
		private final LinkedList<Helix> subCompHelices;
		
		public SubOutput(Helix subMainHelix,LinkedList<Helix> subCompHelices) {
			this.subCompHelices = subCompHelices;
			this.subMainHelix = subMainHelix;
		}
		
		public LinkedList<Helix> getSubCompHelices() {
			return subCompHelices;
		}
		
		public Helix getSubMainHelix() {
			return subMainHelix;
		}
	}
	
	private SubOutput getSubCompHelices(int[] subStruc, NucleotideProbsPrecise nucleotideProbs, double weight, int compLeft, int compRight) {
		IOsideCalculator subIOCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		int[][] subDist = new DistancesCalculator().distCalc(subStruc);
		boolean[][] subCanPair = new PossiblePairFinder().canPair(subStruc);
		InsideOutsideProbabilities subInside = subIOCalc.inside(nucleotideProbs, subDist, weight, subStruc, subCanPair);
		InsideOutsideProbabilities subOutside = subIOCalc.outside(subInside, nucleotideProbs, subDist, weight, subStruc, subCanPair);
		PosteriorProbabilitiesCalculator subPosteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities subPosterior = subPosteriorProbabilitiesCalculator.calculate(subInside, subOutside, 
				nucleotideProbs, subDist, weight, subStruc, subCanPair);
		BigDecimal[][] subDiffs = subPosteriorProbabilitiesCalculator.getDiffs(subPosterior.getPairedProbs(), subPosterior.getUnpairedProbs(), subCanPair);
		Helix subMainHelix = new HelicesMaker().makeHelix(compLeft, compRight, subDiffs, subCanPair); 
		boolean[][] subIncomp = new IncompatiblePairsFinder().findAll(subCanPair, subMainHelix);
		LinkedList<Helix> subCompHelices = new IncompatiblePairsFinder().getIncompatibleHelices(subCanPair, subIncomp, 
				subDiffs, compLeft, compRight);
		SubOutput subOutput = new SubOutput(subMainHelix,subCompHelices);
		return subOutput;
	}
	
	private SubOutput getSubCompHelicesE(int[] subStruc, NucleotideProbsPrecise nucleotideProbs, int compLeft, int compRight) {
		IOsideCalculator subIOCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		boolean[][] subCanPair = new PossiblePairFinder().canPair(subStruc);
		InsideOutsideProbabilities subInside = subIOCalc.insideE(nucleotideProbs, subStruc, subCanPair);
		InsideOutsideProbabilities subOutside = subIOCalc.outsideE(subInside, nucleotideProbs, subStruc, subCanPair);
		PosteriorProbabilitiesCalculator subPosteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities subPosterior = subPosteriorProbabilitiesCalculator.calculateE(subInside, subOutside, 
				nucleotideProbs, subStruc, subCanPair);
		BigDecimal[][] subDiffs = subPosteriorProbabilitiesCalculator.getDiffs(subPosterior.getPairedProbs(), subPosterior.getUnpairedProbs(), subCanPair);
		Helix subMainHelix = new HelicesMaker().makeHelix(compLeft, compRight, subDiffs, subCanPair); 
		boolean[][] subIncomp = new IncompatiblePairsFinder().findAll(subCanPair, subMainHelix);
		LinkedList<Helix> subCompHelices = new IncompatiblePairsFinder().getIncompatibleHelices(subCanPair, subIncomp, 
				subDiffs, compLeft, compRight);
		SubOutput subOutput = new SubOutput(subMainHelix,subCompHelices);
		return subOutput;
	}
	
	public PPOutput calculateE(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateE(insideProbs, outsideProbs, nucleotideProbs, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		BigDecimal[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		BigDecimal[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutput output = null;		
		Helix chosenHelix = null;
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutput(-1,-1,0,BigDecimal.ZERO,BigDecimal.ZERO);
		} else {
			BigDecimal[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			BigDecimal diff = diffs[leftIdx][rightIdx];
			Helix mainHelix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			BigDecimal rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			int mainRight = mainHelix.getRightIdx();
			BigDecimal mainScore = mainHelix.getScore();
			// test whether we should be more careful about our choice of pairing...
			if (rprob.compareTo(Constants.DynamicCutOff)<0) {
				incomp = new IncompatiblePairsFinder().findAll(canPair, mainHelix);
				System.out.println("Possible incompatibility (" + rprob.doubleValue() + ") detected... Proposed new structure is (score " + mainScore.doubleValue() +"):");
				int[] tmpstructure = new StructureUtils().makeNewStructure(structure, mainHelix);
				System.out.println(new LoggingOutputGenerator().dumpStructure(tmpstructure));
				LinkedList<Helix> incompHelices = new IncompatiblePairsFinder().getIncompatibleHelices(canPair,incomp,diffs,leftIdx,rightIdx);
				LinkedList<Helix> interestingHelices = new LinkedList<Helix>();
				System.out.println("Have determined " + incompHelices.size() + " incompatible helices.");
				for (Helix competingHelix: incompHelices) {
					System.out.println("Testing alternative structure: (score " + competingHelix.getScore().doubleValue() +"):");
					tmpstructure = new StructureUtils().makeNewStructure(structure, competingHelix);
					System.out.println(new LoggingOutputGenerator().dumpStructure(tmpstructure));
					int compRight = competingHelix.getRightIdx();
					if (competingHelix.getScore().compareTo(mainScore)>0) { //TODO: replace this with some threshold
						if (mainRight < compRight) {
							System.out.println("Competing structure scores more highly and competing helix further downstream. Breaking proposed helix.");
							// "break" main helix and form competing Helix
							interestingHelices.add(competingHelix);
						} else {
							System.out.println("Competing structure upstream from proposed structure, but scores more highly. Testing for transient structures.");
							// see whether we can prevent the competing Helix from forming
							int compLeft = competingHelix.getLeftIdx();
							int[] subStruc = new StructureUtils().getSubStructure(structure,0,compRight+Constants.TranscriptionOffset); 
							// now get Helices competing with Substructure
							SubOutput subOutput = getSubCompHelicesE(subStruc, nucleotideProbs, compLeft, compRight);
							Helix subMainHelix = subOutput.getSubMainHelix();
							LinkedList<Helix> subCompHelices = subOutput.getSubCompHelices();
							BigDecimal subMainScore = subMainHelix.getScore();
							boolean found = false;
							for (Helix subCompHelix: subCompHelices) {
								if (subCompHelix.getScore().compareTo(subMainScore)>0) {
									found = true; break;
								}
							}
							if (!found) {
								interestingHelices.add(competingHelix);
								System.out.println("No transient structure found. Breaking main helix.");
							} else {
								System.out.println("Found transient structure. Breaking competing helix.");
							}
						}
					} else if (mainRight < compRight) {
						System.out.println("Proposed helix scores more highly than competing helix downstream. Testing for transient structures.");
						// see whether we can prevent the main Helix from forming
						int mainLeft = mainHelix.getLeftIdx();
						int[] subStruc = new StructureUtils().getSubStructure(structure,0,mainRight+Constants.TranscriptionOffset); 
						// now get matrix of probabilities for substructure
						SubOutput subOutput = getSubCompHelicesE(subStruc, nucleotideProbs, mainLeft, mainRight);
						Helix subMainHelix = subOutput.getSubMainHelix();
						LinkedList<Helix> subCompHelices = subOutput.getSubCompHelices();
						BigDecimal subMainScore = subMainHelix.getScore();
						boolean found = false;
						for (Helix subCompHelix: subCompHelices) {
							if (subCompHelix.getScore().compareTo(subMainScore)>0) {
								found = true; break;
							}
						}
						if (found) {
							interestingHelices.add(competingHelix);
							System.out.println("Transient structure found. Breaking proposed helix.");
						} else {
							System.out.println("No transient structure found. Sticking with proposed helix.");
						}
						// otherwise "break" competing helix, i.e. do not include competing helix in list of interesting helices
					}
				}
				// if no interesting structures: stick with main structure; otherwise compare interesting competing structures and decide.
				if (!interestingHelices.isEmpty()) {
					BigDecimal maxScore = BigDecimal.ZERO;
					chosenHelix = mainHelix;
					for (Helix helix: interestingHelices) {
						if (maxScore.compareTo(helix.getScore())<0) {
							chosenHelix = helix;
							maxScore = helix.getScore();
						}
					}	
				} else {
					chosenHelix = mainHelix;
				}
			} else {
				chosenHelix = mainHelix;
			}
			output = new PPOutput(chosenHelix.getLeftIdx(), chosenHelix.getRightIdx(), chosenHelix.getHelixLength(), diff, rprob);			
		}
				
		return output;
	}
}
