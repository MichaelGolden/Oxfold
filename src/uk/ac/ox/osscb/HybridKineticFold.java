package uk.ac.ox.osscb;

import java.util.HashMap;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResultDouble;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutputDouble;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ShortDoublePosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class HybridKineticFold {

	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		if(weight <= 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Zero weight not supported at the moment. Input: %f", weight));
		}
		
		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		
		boolean[] blockedColumns = getBlockedColumns(align,parameters.getSAlphabet().getSynonyms());
		
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
		boolean exitBecauseOfDiff = false;
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		ShortDoublePosteriorProbabilities currentPostProbs = ppCalc.calculateShortE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
					
			PPOutputInternalResultDouble postProbs = evol ? 
					kFPppCalc.calculateHybridPpOutput(alignmentProbsEvol, structure,alignmentFile,currentPostProbs,blockedColumns) 
						:kFPppCalc.calculateHybridPpOutput(alignmentProbs, structure,alignmentFile,currentPostProbs,blockedColumns);
					
			PPOutputDouble ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();

			
			//if (ppProbs.getDiff().signum()>0) {
			if ((!evol)&&(ppProbs.getDiff()>Constants.EndNonEvolutionaryFoldDouble)) {
				evol = true;
				continue;
			}
			if (ppProbs.getDiff()>Constants.IterationCutOffDouble) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs.getHelix());
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
	
	
	private boolean[] getBlockedColumns(String[] align, HashMap<String,String[]> synonyms) {
		int length = align[0].length();
		int seqnum = align.length;
		boolean[] blockedColumns = new boolean[length];
		for (int i = 0; i<length;i++) {
			int unsure = 0;
			for (int j = 0; j<seqnum; j++) {
				String obs = align[j].substring(i, i+1);
				if (synonyms.get(obs).length>1) {
					unsure++;
				}
			}
			if (unsure>=0.5*seqnum) {
				blockedColumns[i] = true;
			}
		}
		return blockedColumns;
	}
	
	private void dumpExitReason(boolean exitBecauseOfDiff) {
		ProgramOutput.outMsg(String.format("\texiting because of %s", 
				exitBecauseOfDiff ? "difference under threshold." : "number of iterations."));
	}

	private void dumpCurrentOutput(PPOutputDouble ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g"
				, ppProbs.getLeftIdx()
				, ppProbs.getRightIdx()
				, ppProbs.gethelixLength()
				, ppProbs.getDiff()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}

	
}
