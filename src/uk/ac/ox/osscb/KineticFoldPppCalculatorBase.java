package uk.ac.ox.osscb;

import org.slf4j.Logger;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPOutputDouble;
import uk.ac.ox.osscb.inoutside.PPOutputHelix;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.ShortDoublePosteriorProbabilities;

/**
 * Base for both {@link KineticFoldPppCalculatorWeightLess} and
 * {@link KineticFoldPppCalculatorWithWeight} in order to log inside/outside
 * matrices upon <b>first</b> iteration only. 
 * Not sure it's a best design  (base class knows about some details 
 * about how subclasses work - inside/outside) 
 * 
 * @author Vladimir
 */
public abstract class KineticFoldPppCalculatorBase implements KineticFoldPppCalculator {

	protected final Grammar grammar;
	protected final IOsideCalculator ioCalc;
	private int iterNo = 0;
	
	public KineticFoldPppCalculatorBase(Grammar grammar, IOsideCalculator ioCalc) {
		super();
		this.grammar = grammar;
		this.ioCalc = ioCalc;
	}
	
	/**
	 * we use implementation-specific logger. 
	 * @return
	 */
	protected abstract Logger getLog();


	public PPOutput calculatePpOutput(NucleotideProbsPrecise alignmentProbs,
			int[] structure) {
		
		this.iterNo++;
		//still returns PPOutput as required
		PPOutputInternalResult2 calculatePpOutputInternal = calculatePpOutputInternal(alignmentProbs, structure, null);
		
		if(1 == this.iterNo){
			// print the matrices
			// log.info(String.format());
		}
		
		return calculatePpOutputInternal.getPpProbs();
	}
	
	@Override
	public PPOutputInternalResult2 calculatePpOutputInternalResult2(NucleotideProbsPrecise alignmentProbs,
			int[] structure, PosteriorProbabilities postProbs) {
		
		this.iterNo++;
		
		PPOutputInternalResult2 calculatePpOutputInternal = calculatePpOutputInternal(alignmentProbs, structure, postProbs);
		
		if(1 == this.iterNo){
			// print the matrices
			// log.info(String.format());
		}
		
		return calculatePpOutputInternal;
	}
	
	public abstract PPOutputInternalResult2 calculatePpOutputInternal(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, PosteriorProbabilities postProbs);

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

