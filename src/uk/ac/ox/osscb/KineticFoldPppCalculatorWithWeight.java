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
 * high-level {@link PPOutput} calculation. Weight <b>is</b> used
 * 
 * @author Vladimir
 *
 */
public class KineticFoldPppCalculatorWithWeight 
	extends KineticFoldPppCalculatorBase {
	
	private final Logger log = LoggerFactory.getLogger(KineticFoldPppCalculatorWithWeight .class);
	
	private final double weight;

	public KineticFoldPppCalculatorWithWeight(double weight, Grammar grammar,
			IOsideCalculator ioCalc) {
		super(grammar, ioCalc);
		this.weight = weight;
	}
	
	@Override
	protected Logger getLog() {
		return this.log;
	}
	
	public PPOutputInternalResult2 calculatePpOutputInternal(NucleotideProbsPrecise alignmentProbs, int[] structure, PosteriorProbabilities prevProbs) {
		
		double[][] distances = new DistancesCalculator2().distCalc(structure, prevProbs);
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, weight, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, weight, structure, canPair);				
		PPOutput ppProbs = new PPProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		//PosteriorProbabilities completePPProbs = new PosteriorProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		PosteriorProbabilities completePPProbs = new PosteriorProbabilitiesCalculator(grammar).calculateParallel(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);

		
		return new PPOutputInternalResult2(ppProbs, completePPProbs);
	}
	
	public PPOutput calculateDynamicPpOutput(NucleotideProbsPrecise alignmentProbs, int[] structure, boolean[][] canPair) {
		int[][] distances = new DistancesCalculator().distCalc(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, weight, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, weight, structure, canPair);				
		PPOutput ppProbs = new DynamicPPProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		return ppProbs;
	}
	
	public PPOutputHelixInternalResult calculateCotranscriptionalPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair, PosteriorProbabilities posteriorProbs) {
		
		double[][] distances = new DistancesCalculator2().distCalc(structure, posteriorProbs);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, weight, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, weight, structure, canPair);				
		PPOutputHelix ppProbs = new PPProbabilitiesCalculator(grammar).calculateCT(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		PosteriorProbabilities completePPProbs = new PosteriorProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		
		return new PPOutputHelixInternalResult(ppProbs, completePPProbs);
	}
	
	public PPOutputInternalResultDouble calculateHybridPpOutput(NucleotideProbsPrecise alignmentProbs, int[] structure, String alignmentPath, 
			ShortDoublePosteriorProbabilities postProbs, boolean[] blockedColumns) {
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		double[][] distances = new DistancesCalculator2().distCalc(structure, postProbs);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, weight, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, weight, structure, canPair);				
		PPOutputDouble ppProbs = new PPProbabilitiesCalculator(grammar).calculateHybrid(insideProbs,outsideProbs,alignmentProbs,distances,weight,structure,canPair,alignmentPath,blockedColumns);
		ShortDoublePosteriorProbabilities completePPProbs = new PosteriorProbabilitiesCalculator(grammar).calculateShort(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
		return new PPOutputInternalResultDouble(ppProbs, completePPProbs);
	}



	
}
