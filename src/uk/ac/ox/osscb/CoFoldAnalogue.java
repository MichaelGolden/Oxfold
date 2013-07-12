package uk.ac.ox.osscb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputHelixInternalResult;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.CoFoldInsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.CoFoldPosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPOutputHelix;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class CoFoldAnalogue {
	private final Logger log = LoggerFactory.getLogger(CoFoldAnalogue.class);
	
	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double alpha, double tau){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		/*if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}*/

		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);
		
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
		
		CoFoldInsideOutsideCalculator ioCalc = new CoFoldInsideOutsideCalculator(grammar);
		
		
		//IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new CoFoldInsideOutsideCalculator(grammar), grammar);
		//alpha = 0.5;
		//tau = 640;
		CoFoldPppCalculator cofoldCalc = new CoFoldPppCalculator(alpha, tau,grammar,ioCalc);
		/*KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);*/
		
		//ioCalc.outside(insideProbs, pairingProbs, distances, alpha, tau, structure, canPair)
		
		boolean evol = true;
		boolean exitBecauseOfDiff = false;
	
		boolean [][] canPair = new boolean [structure.length][structure.length];
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = 0 ; j < canPair.length ; j++)
			{
				canPair[i][j] = true;
			}
		}
		InsideOutsideProbabilities insideProbs = ioCalc.inside(evol ? alignmentProbs : alignmentProbsEvol, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, alpha, tau, structure, canPair);
		CoFoldPosteriorProbabilitiesCalculator ppCalc = new CoFoldPosteriorProbabilitiesCalculator(grammar);
		int [][] distances = null;
		PosteriorProbabilities currentPostProbs = ppCalc.calculate(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		MEACalculator meaCalculator = new MEACalculator();
		structure = meaCalculator.calculate(currentPostProbs);
		//System.out.println(RNAFoldingTools.);
		
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		//InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		//InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		//PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		//PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		/*
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			canPair = new PossiblePairFinder().canPair(structure);			
					
			PPOutputInternalResult2 postProbs = evol ? 
					cofoldCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						: cofoldCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();

			//if (ppProbs.getDiff().signum()>0) {
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
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			dumpStructure(structure);
		}*/


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName(), structure);
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
				, ppProbs.getComp()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}
}
