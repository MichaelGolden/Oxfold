package uk.ac.ox.osscb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputHelixInternalResult;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResultDouble;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.inoutside.CoFoldInsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.CoFoldPosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.CoFoldProbabilitiesCalculator;
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
public class CoFoldPppCalculator {
	//extends KineticFoldPppCalculatorBase {
	
	private final Logger log = LoggerFactory.getLogger(CoFoldPppCalculator.class);
	
	private final double alpha;
	private final double tau;
	
	Grammar grammar;
	CoFoldInsideOutsideCalculator ioCalc;
	
	public CoFoldPppCalculator(double alpha, double tau, Grammar grammar,
			CoFoldInsideOutsideCalculator ioCalc) {
		this.grammar = grammar;
		this.ioCalc = ioCalc;
		//super(grammar, ioCalc);
		this.alpha = alpha;
		this.tau = tau;
	}
	
	//@Override
	protected Logger getLog() {
		return this.log;
	}
	
	
	public PPOutputInternalResult2 calculatePpOutputInternal(NucleotideProbsPrecise alignmentProbs, int[] structure, PosteriorProbabilities prevProbs) {
		
		double[][] distances = new DistancesCalculator2().distCalc(structure, prevProbs);
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);				
		PPOutput ppProbs = new CoFoldProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		PosteriorProbabilities completePPProbs = new CoFoldPosteriorProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);

		
		return new PPOutputInternalResult2(ppProbs, completePPProbs);
	}
	
	public PPOutput calculateDynamicPpOutput(NucleotideProbsPrecise alignmentProbs, int[] structure, boolean[][] canPair) {
		int[][] distances = new DistancesCalculator().distCalc(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);				
		PPOutput ppProbs = new CoFoldProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		return ppProbs;
	}
	
	public PPOutputHelixInternalResult calculateCotranscriptionalPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair, PosteriorProbabilities posteriorProbs) {
		
		double[][] distances = new DistancesCalculator2().distCalc(structure, posteriorProbs);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);				
		PPOutputHelix ppProbs = new PPProbabilitiesCalculator(grammar).calculateCT(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		PosteriorProbabilities completePPProbs = new CoFoldPosteriorProbabilitiesCalculator(grammar).calculate(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		
		return new PPOutputHelixInternalResult(ppProbs, completePPProbs);
	}
	
	public PPOutputInternalResultDouble calculateHybridPpOutput(NucleotideProbsPrecise alignmentProbs, int[] structure, String alignmentPath, 
			ShortDoublePosteriorProbabilities postProbs, boolean[] blockedColumns) {
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		double[][] distances = new DistancesCalculator2().distCalc(structure, postProbs);
		InsideOutsideProbabilities insideProbs = ioCalc.inside(alignmentProbs, distances, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);				
		PPOutputDouble ppProbs = new CoFoldPppCalculator(grammar).calculateHybrid(insideProbs,outsideProbs,alignmentProbs,distances,alpha, tau,structure,canPair,alignmentPath,blockedColumns);
		ShortDoublePosteriorProbabilities completePPProbs = new PosteriorProbabilitiesCalculator(grammar).calculateShort(insideProbs, outsideProbs, alignmentProbs, distances, alpha, tau, structure, canPair);
		return new PPOutputInternalResultDouble(ppProbs, completePPProbs);
	}

	protected class PPOutputInternalResult{
		private InsideOutsideProbabilities insideProbs;
		private InsideOutsideProbabilities outsideProbs;
		private PPOutput ppProbs;
		public PPOutputInternalResult(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, PPOutput ppProbs) {
			super();
			this.insideProbs = insideProbs;
			this.outsideProbs = outsideProbs;
			this.ppProbs = ppProbs;
		}
		public InsideOutsideProbabilities getInsideProbs() {
			return insideProbs;
		}
		public InsideOutsideProbabilities getOutsideProbs() {
			return outsideProbs;
		}
		public PPOutput getPpProbs() {
			return ppProbs;
		}
	}

	protected class PPOutputInternalResult2{
		private PPOutput ppProbs;
		private PosteriorProbabilities posteriorProbs;
		public PPOutputInternalResult2(PPOutput ppProbs,
				PosteriorProbabilities posteriorProbs) {
			super();
			this.ppProbs = ppProbs;
			this.posteriorProbs = posteriorProbs;
		}
		
		public PPOutput getPpProbs() {
			return ppProbs;
		}
		
		public PosteriorProbabilities getPosteriorProbs(){
			return posteriorProbs;
		}
	}
	
	protected class PPOutputInternalResultDouble{
		private PPOutputDouble ppProbs;
		private ShortDoublePosteriorProbabilities posteriorProbs;
		public PPOutputInternalResultDouble(PPOutputDouble ppProbs,
				ShortDoublePosteriorProbabilities posteriorProbs) {
			super();
			this.ppProbs = ppProbs;
			this.posteriorProbs = posteriorProbs;
		}
		public PPOutputDouble getPpProbs() {
			return ppProbs;
		}
		public ShortDoublePosteriorProbabilities getPosteriorProbs(){
			return posteriorProbs;
		}
	}
	
	protected class PPOutputHelixInternalResult{
		private PPOutputHelix ppProbs;
		private PosteriorProbabilities posteriorProbs;
		public PPOutputHelixInternalResult(PPOutputHelix ppProbs,
				PosteriorProbabilities posteriorProbs) {
			super();
			this.ppProbs = ppProbs;
			this.posteriorProbs = posteriorProbs;
		}
		public PPOutputHelix getPpProbs() {
			return ppProbs;
		}
		public PosteriorProbabilities getPosteriorProbs(){
			return posteriorProbs;
		}
	}
	

	
}
