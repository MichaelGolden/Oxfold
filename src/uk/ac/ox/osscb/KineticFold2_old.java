package uk.ac.ox.osscb;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

/**
 * Functional entry point. Does not parse input arguments, etc.
 * 
 * @author Vladimir, lepuslapis
 *
 */
public class KineticFold2_old {
	
	private final Logger log = LoggerFactory.getLogger(KineticFold2_old.class);


	/**
	 * 
	 * @param alignmentFile
	 * @param grammarFile
	 * @param paramsFile
	 * @param weight
	 */
	public void fold(String alignmentFile, String grammarFile, String paramsFile, double weight){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}
		
		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParser2().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile, parameters.getSAlphabet());
				
		// int[][] alignment = new AlignmentConverter().convert(align, parameters);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters); 

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
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities completePPProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilities originalProbs = completePPProbs;
		double[][] distances = new DistancesCalculator2().distCalc(structure, completePPProbs);
		PosteriorProbabilities currentPostProbs = completePPProbs;
		
		boolean exitBecauseOfDiff = false;
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			canPair = new PossiblePairFinder().canPair(structure);				
			completePPProbs = ppCalc.calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
			
			PPOutputInternalResult2 postProbs = kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();
			
			//if (ppProbs.getDiff().signum()>0) {
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			OutputGenerator outputGenerator = new LoggingOutputGenerator();
			outputGenerator.generate(structure);
			dumpStructure(structure);
		}
		

		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			
			probs.savePosteriorProbabilities(new File(alignmentFile+".noevol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		writeDotBracketFile(new File(alignmentFile+".noevol.dbn"),new File(alignmentFile).getName(), structure);
	}
	
	public static void writeDotBracketFile(File dbnFile, String title, int [] structure)
	{
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(dbnFile));
			writer.write(">"+title+"\n");
			writer.write("\n");
			writer.write(LoggingOutputGenerator.dumpStructure(structure)+"\n");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight, double delta2){
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
		
		boolean [] delete =new boolean[align[0].length()];
		//boolean [] delete = CoFoldAnalogue.getGappyColumns(align, 0.25);
		align = CoFoldAnalogue.deleteColumns(align, delete);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);
		
		// Added in
		//
		
		
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		alignmentProbsEvol = EvolProbabilitiesCalculator.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile));
		
		
		if(log.isDebugEnabled()){
			log.debug(String.format("Pairing probs (NON-evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (NON-evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
			
			log.debug(String.format("Pairing probs (Evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (Evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
		}
		
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
		boolean exitBecauseOfDiff = false;
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure, delete);
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		int fastIterations = 0;
		int iterSoFar = 0;
		for(; iterSoFar < Constants.MaxIterations; iterSoFar++){
			//canPair = new PossiblePairFinder().canPair(structure);		
			canPair = new PossiblePairFinder().canPair(structure, delete);
					
			PPOutputInternalResult2 postProbs = evol ? 
					kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();
			
			/*try {
				currentPostProbs.savePosteriorProbabilities(new File("iteration"+iterSoFar+".currentPostProbs"));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/

			//if (ppProbs.getDiff().signum()>0) {
			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				
				double currentDelta = ppProbs.getDiff().doubleValue();
				boolean runfast = true;
				for(int c = 0 ; runfast ; c++)
				{
					fastIterations++;
					if(ppProbs.getDiff().compareTo(new PointRes(delta2))<0)
					{
						break;
					}
					boolean [][] canPair2 = new PossiblePairFinder().canPair(structure);					

					PointRes[] unpairedProbs = currentPostProbs.getUnpairedProbs();
					PointRes[][] pairedProbs = currentPostProbs.getPairedProbs();
					
					int leftIdx = -1;
					int rightIdx = -1;
					
					for(int i = 0 ; i < canPair2.length ; i++)
					{
						for(int j = i + 1 ; j < canPair2.length ; j++)
						{
							if(canPair2[i][j]  && structure[i] < 0 && structure[j] < 0)
							{
								if((leftIdx == -1 && rightIdx == -1) || pairedProbs[i][j].compareTo( pairedProbs[leftIdx][rightIdx]) > 0)
								{
									leftIdx = i;
									rightIdx = j;
								}
							}
						}
					}
					System.out.println("fast "+c+"\t"+leftIdx+"\t"+rightIdx);
					
					if ((leftIdx<0)||(rightIdx<0)) {
						ppProbs = new PPOutput(-1,-1,0,PointRes.ZERO,PointRes.ZERO);
					} else {
						boolean[][] incomp = new IncompatiblePairsFinder().find(canPair2, leftIdx, rightIdx);
						PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
						PointRes[][] diffs = PosteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair2);
						PointRes diff = diffs[leftIdx][rightIdx];
						Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair2);
						PPOutput ppProbs2 = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
	
						//if(ppProbs2.getDiff().compareTo(new PointRes(Math.max(delta2,weight)))>0)
						if(ppProbs2.getDiff().compareTo(new PointRes(Math.max(currentDelta*delta2,weight)))>0)
						{
							canPair = canPair2;
							ppProbs = ppProbs2;
							int [] structure2 = new int[structure.length];
							for(int i = 0 ; i < structure.length ; i++)
							{
								structure2[i] = structure[i]+1;
							}
							//System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(structure2));
							structure = new StructureUtils().makeNewStructure(structure, ppProbs2);
	
							continue;
						}
					}
					break;				
				}
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			dumpStructure(structure);
		}


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			// Added in
			probs.pairedProbs = CoFoldAnalogue.reinsertDeleted(currentPostProbs.pairedProbs, delete, PointRes.ZERO);
			structure = CoFoldAnalogue.reinsertDeleted(structure, delete, Constants.UnpairedBaseIdx);
			probs.savePosteriorProbabilities(new File(alignmentFile+".evol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String metadata = ";oxfolditer="+iterSoFar+",fastiter="+fastIterations+",delta2="+delta2;
		System.out.println(metadata);
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName()+metadata, structure);
	}
	
	private void dumpExitReason(boolean exitBecauseOfDiff) {
		ProgramOutput.outMsg(String.format("\texiting because of %s", 
				exitBecauseOfDiff ? "difference under threshold." : "number of iterations."));
	}

	private void dumpCurrentOutput(PPOutput ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g\trprob: %g"
				, ppProbs.getLeftIdx()
				, ppProbs.getRightIdx()
				, ppProbs.gethelixLength()
				, ppProbs.getDiff().doubleValue()
				, ppProbs.getComp().doubleValue()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}
}
