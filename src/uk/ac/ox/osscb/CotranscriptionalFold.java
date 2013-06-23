package uk.ac.ox.osscb;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.LinkedList;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputHelixInternalResult;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPOutputHelix;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class CotranscriptionalFold {
	
	/**
	 * Co-transcriptional Folding
	 * @param alignmentFile
	 * @param grammarFile
	 * @param paramsFile
	 * @param treeFile
	 * @param weight
	 * @author Lepuslapis
	 */
	
	public void newFoldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight) {
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		if (weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}

		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);
				
		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		
		IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);

		int seqEnd = Constants.TranscriptionJump-1;
		int[] transStruc = new StructureUtils().getSubStructure(structure, 0, seqEnd);
		boolean transcriptionEnd = false;
		
		for (int transIter = 1; transIter < Constants.MaxIterations; transIter++) {
			boolean evol = false;
			//init distance etc. use transStruc for structure
			transStruc = new StructureUtils().getSubStructure(structure, 0, seqEnd);

			boolean[][] canPair = new PossiblePairFinder().canPair(transStruc);
			InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, transStruc, canPair);
			InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, transStruc, canPair);
			PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
			PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, transStruc, canPair);
			
			for (int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++) {
				canPair = new PossiblePairFinder().canPair(transStruc);
				PPOutputHelixInternalResult postProbs = evol ? kFPppCalc.calculateCotranscriptionalPpOutput(alignmentProbsEvol, transStruc, canPair, currentPostProbs) 
						: kFPppCalc.calculateCotranscriptionalPpOutput(alignmentProbs, transStruc, canPair, currentPostProbs);
				PPOutputHelix ppProbs = postProbs.getPpProbs();
				currentPostProbs = postProbs.getPosteriorProbs();
				
				if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
					evol = true;
					continue;
				}
				if ((ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0)) { //&&(ppProbs.getHelix().getHelixLength()>=Constants.TranscriptionMinimumHelixLength)) {
					transStruc = new StructureUtils().makeNewStructure(transStruc, ppProbs.getHelix());
					updateProbs(alignmentProbsEvol,transStruc.length,ppProbs.getHelix()); updateProbs(alignmentProbs,transStruc.length,ppProbs.getHelix());
					dumpCurrentOutput(ppProbs);
				} else {
					dumpCurrentOutput(ppProbs);
					break;
				}
				OutputGenerator outputGenerator = new LoggingOutputGenerator();
				outputGenerator.generate(transStruc);
				dumpStructure(transStruc);
			}
			System.out.println("Finished folding transcribed structure.");
			//dumpStructure(transStruc);
			if (transcriptionEnd) {
				break;
			}
			seqEnd += Constants.TranscriptionJump;
			if (seqEnd>=structure.length) {
				seqEnd = structure.length-1;
				transcriptionEnd = true;
			}
		}
	}
	
	private void updateProbs(NucleotideProbsPrecise alignmentProbs, int length, Helix helix) {
		int helixLength = helix.getHelixLength();
		int left = helix.getLeftIdx(); int right = helix.getRightIdx();
		double scale = (double) length*length/helixLength;
		BigDecimal priorWeight = BigDecimal.valueOf(18*scale);
		BigDecimal priorWeightMinusOne = priorWeight.subtract(BigDecimal.ONE);
		for (int j = 0; j<helixLength; j++) {
			for (int k = 0; k<left+j; k++) {
				BigDecimal prob1 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(k,left+j)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(k, left+j, prob1);
				BigDecimal prob2 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(k,right-j)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(k, right-j, prob2);
			}
			for (int k = left+j+1; k<right-j; k++) {
				BigDecimal prob1 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(left+j,k)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(left+j, k, prob1);
				BigDecimal prob2 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(k,right-j)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(k,right-j, prob2);
			}
			for (int k = right-j+1;k<length; k++) {
				BigDecimal prob1 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(left+j,k)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(left+j, k, prob1);
				BigDecimal prob2 = priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(right-j,k)).divide(priorWeight,RoundingMode.HALF_UP);
				alignmentProbs.setPairingProbability(right-j, k, prob2);
			}
			BigDecimal prob = (priorWeightMinusOne.multiply(alignmentProbs.getPairingProbability(left+j,right-j)).add(BigDecimal.ONE)).divide(priorWeight,RoundingMode.HALF_UP);
			alignmentProbs.setPairingProbability(left+j, right-j, prob);
		}
	}
	
	/*
	 * Older version of Cotranscriptional Folding Algorithm follows below.
	 */
	/*
	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}

		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);
		
		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		
		//IterationsGenerator iterationsGenerator = new IterationsGenerator ();
		
		IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);
		
		boolean evol = false;
		// start transcription
		int seqEnd = Constants.TranscriptionJump-1;
		int[] transStruc = new int[]{-1};
		LinkedList<Helix> helixList = new LinkedList<Helix>();
		while (seqEnd<structure.length) {
			evol = false;
			transStruc = new StructureUtils().getSubStructure(structure, 0, seqEnd);
			boolean[][] canPair = new PossiblePairFinder().canPair(transStruc);
			for (int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++) {
				PPOutputHelix ppProbs = evol ? kFPppCalc.calculateCotranscriptionalPpOutput(alignmentProbsEvol, transStruc, canPair) 
						: kFPppCalc.calculateCotranscriptionalPpOutput(alignmentProbs, transStruc, canPair);
				// update scores of existing helices
				updateScores(helixList,ppProbs.getDiffs());
				if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
					evol = true;
					continue;
				}
				if ((ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0)&&(ppProbs.getHelix().getHelixLength()>=Constants.TranscriptionMinimumHelixLength)) {
					Helix propHelix = ppProbs.getHelix();
					if (!helixList.contains(propHelix)) {
						BigDecimal oppScore = BigDecimal.ZERO;
						LinkedList<Helix> oppHelices = new LinkedList<Helix>();
						for (Helix oldHelix: helixList) {
							if (interfere(propHelix,oldHelix)) {
								oppScore = oppScore.add(oldHelix.getScore());
								oppHelices.add(oldHelix);
							}
						}
						System.out.println("Have found " + oppHelices.size() + " competing helices among " + helixList.size() + " existing helices.");
						if ((propHelix.getScore().subtract(oppScore)).compareTo(Constants.BreakingCutOff)>0) {
							System.out.println("New Helix beats competing helices.");
							for (Helix helix: oppHelices) {
								helixList.remove(helix);
							}
							helixList.add(propHelix);
							transStruc = new StructureUtils().makeNewStructure(transStruc,propHelix);
						} else {
							System.out.println("Sticking with existing helices.");
							break;
						}
					} else {
						transStruc = new StructureUtils().makeNewStructure(transStruc, propHelix);
					}
					dumpCurrentOutput(ppProbs);
				} else {
					dumpCurrentOutput(ppProbs);
					break;
				}
				dumpStructure(transStruc);
			}
			for (Helix helix: helixList) {
				transStruc = new StructureUtils().makeNewStructure(transStruc, helix);
			}
			System.out.println("Finished folding transcribed structure. More transcription follows.");
			dumpStructure(transStruc);
			seqEnd += Constants.TranscriptionJump;
		}
		updateStructure(structure, transStruc);
		
		evol = false;
		boolean exitBecauseOfDiff = false;
		// Transcription has ended. Now add in shorter helices. 
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			
			PPOutput ppProbs = evol ? kFPppCalc.calculatePpOutput(alignmentProbsEvol, structure) 
					: kFPppCalc.calculatePpOutput(alignmentProbs, structure);
			
			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			dumpStructure(structure);
		}


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);

	} */
	
	private boolean interfere(Helix helixA, Helix helixB) {
		int leftA = helixA.getLeftIdx(); int leftB = helixB.getLeftIdx();
		int rightA = helixA.getRightIdx(); int rightB = helixB.getRightIdx();
		int lengthA = helixA.getHelixLength(); int lengthB = helixB.getHelixLength();
		if (rightA<leftB) {
			return false;
		} else if (leftA>rightB) {
			return false;
		} else if ((leftA+lengthA<=leftB)&&(rightA-lengthA>=rightB)) {
			return false;
		} else if ((leftB+lengthB<=leftA)&&(rightB-lengthB>=rightA)) {
			return false;
		} else {
			return true;
		}		
	}
	
	private void updateScores(LinkedList<Helix> helixList, BigDecimal[][] diffs) {
		for (Helix helix: helixList) {
			BigDecimal scoretmp = BigDecimal.ZERO;
			int left = helix.getLeftIdx(); int right = helix.getRightIdx();
			for (int j = 0; j<helix.getHelixLength(); j++) {
				scoretmp = scoretmp.add(diffs[left+j][right-j]);
			}
			helix.setScore(scoretmp);
		}
	}
	
	private void updateStructure(int[] structure, int[] transStruc) {
		for (int j = 0; j<transStruc.length; j++) {
			structure[j] = transStruc[j];
		}
	}
	
	
	
	private void dumpExitReason(boolean exitBecauseOfDiff) {
		ProgramOutput.outMsg(String.format("\texiting because of %s", 
				exitBecauseOfDiff ? "difference under threshold." : "number of iterations."));
	}

	private void dumpCurrentOutput(PPOutputHelix ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g\trprob: %g"
				, ppProbs.getHelix().getLeftIdx()
				, ppProbs.getHelix().getRightIdx()
				, ppProbs.getHelix().getHelixLength()
				, ppProbs.getDiff().doubleValue()
				, ppProbs.getComp()
				);
		ProgramOutput.outMsg(msg);
	}
	
	private void dumpCurrentOutput(PPOutput ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g\trprob: %g"
				, ppProbs.getLeftIdx()
				, ppProbs.getRightIdx()
				, ppProbs.gethelixLength()
				, ppProbs.getDiff().doubleValue()
				, ppProbs.getComp()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}

}
