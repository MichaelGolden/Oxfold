package uk.ac.ox.osscb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.inoutside.DynamicPPProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPOutputDouble;
import uk.ac.ox.osscb.inoutside.PPOutputHelix;
import uk.ac.ox.osscb.inoutside.PPProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ShortDoublePosteriorProbabilities;

/**
 * high-level {@link PPOutput} calculation. Weight is not used
 * 
 * @author Vladimir
 */
public class KineticFoldPppCalculatorWeightLess 
	extends KineticFoldPppCalculatorBase {
	
	private final Logger log = LoggerFactory.getLogger(KineticFoldPppCalculatorWeightLess.class);
	

	public KineticFoldPppCalculatorWeightLess(Grammar grammar,
			IOsideCalculator ioCalc) {
		super(grammar, ioCalc);
	}


	@Override
	public PPOutputInternalResult2 calculatePpOutputInternal(
			NucleotideProbsPrecise alignmentProbs, int[] structure,
			PosteriorProbabilities postProbs) {
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);		
		PosteriorProbabilities posteriorProbs = new PosteriorProbabilitiesCalculator(grammar).calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		PPOutput ppProbs = new PPProbabilitiesCalculator(grammar).calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		return new PPOutputInternalResult2(ppProbs, posteriorProbs);
	}
	
	
	public PPOutput calculateDynamicPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair) {
		
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);				
		PPOutput ppProbs = new DynamicPPProbabilitiesCalculator(grammar).calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		return ppProbs;
	}

	public PPOutputHelixInternalResult calculateCotranscriptionalPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair, PosteriorProbabilities postProbs) {
		
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);				
		PPOutputHelix ppProbs = new PPProbabilitiesCalculator(grammar).calculateCTE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilities posteriorProbs = new PosteriorProbabilitiesCalculator(grammar).calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		return new PPOutputHelixInternalResult(ppProbs, posteriorProbs);
	}
	
	public PPOutputInternalResultDouble calculateHybridPpOutput(NucleotideProbsPrecise alignmentProbs, int[] structure, String alignmentPath, ShortDoublePosteriorProbabilities postProbs,boolean[] blockedColumns) {
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PPOutputDouble ppProbs = new PPProbabilitiesCalculator(grammar).calculateHybridE(insideProbs,outsideProbs,alignmentProbs,structure,canPair,alignmentPath,blockedColumns);
		ShortDoublePosteriorProbabilities posteriorProbs = new PosteriorProbabilitiesCalculator(grammar).calculateShortE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		return new PPOutputInternalResultDouble(ppProbs, posteriorProbs);
	}
	
	@Override
	protected Logger getLog() {
		return this.log;
	}



}
