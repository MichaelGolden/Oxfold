package uk.ac.ox.osscb;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

/**
 * Functional entry point. Does not parse input arguments, etc.
 * 
 * @author Vladimir, lepuslapis
 *
 */
public class KineticFold {
	
	private final Logger log = LoggerFactory.getLogger(KineticFold.class);


	/**
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
		
		boolean exitBecauseOfDiff = false;
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
						
			PPOutput ppProbs = kFPppCalc.calculatePpOutput(alignmentProbs, structure);
			
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
	}

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
		
		if(log.isDebugEnabled()){
			log.debug(String.format("Pairing probs (NON-evol):%n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (NON-evol):%n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
			
			log.debug(String.format("Pairing probs (Evol):%n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (Evol):%n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
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
		
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			
			PPOutput ppProbs = evol ? 
					kFPppCalc.calculatePpOutput(alignmentProbsEvol, structure) 
						:kFPppCalc.calculatePpOutput(alignmentProbs, structure);
			
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
		}


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);

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
