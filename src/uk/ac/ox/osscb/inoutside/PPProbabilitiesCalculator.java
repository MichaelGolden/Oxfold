package uk.ac.ox.osscb.inoutside;


import java.util.InputMismatchException;
import java.util.Locale;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.HelicesMaker;
import uk.ac.ox.osscb.IncompatiblePairsFinder;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.Util;
import uk.ac.ox.osscb.config.Settings;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.vienna.Fasta2StockholmConvertingViennaRnaRunner;

public class PPProbabilitiesCalculator {
	
	private static final Logger log = LoggerFactory.getLogger(PPProbabilitiesCalculator.class);

	private Grammar grammar;
	
	static {
		// output only once. do not use logger config because we want this to be output only once
		// per program invocation
		log.info(String.format("%s %s:",
				Util.dateTime(), PPProbabilitiesCalculator.class.getName()));
	}
	
	public PPProbabilitiesCalculator(Grammar grammar) {
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
			double[][] distances, double weight, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		//PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculate(insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateParallel(insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		PointRes[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		PointRes[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutput output = null;
		
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutput(-1,-1,0,PointRes.ZERO,PointRes.ZERO);
		} else {
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			PointRes[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			PointRes diff = diffs[leftIdx][rightIdx];
			//Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,pairedProbs,diffs,canPair);
			output = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
			
			logResults(leftIdx, rightIdx, diff, rprob);			
		}

		return output;
	}
	
	public PPOutput calculate(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs,
			int[][] distances, double weight, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculate(insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		PointRes[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		PointRes[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutput output = null;
		
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutput(-1,-1,0,PointRes.ZERO,PointRes.ZERO);
		} else {
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			PointRes[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			PointRes diff = diffs[leftIdx][rightIdx];
			//Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,pairedProbs,diffs,canPair);
			output = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
			
			logResults(leftIdx, rightIdx, diff, rprob);			
		}

		return output;
	}
	
	
	
	public PPOutput calculateE(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateE(insideProbs, outsideProbs, nucleotideProbs, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		PointRes[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		PointRes[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutput output = null;
						
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutput(-1,-1,0,PointRes.ZERO,PointRes.ZERO);
		} else {
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			PointRes[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			PointRes diff = diffs[leftIdx][rightIdx];
			//Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx, pairedProbs,diffs,canPair);
			output = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
			
			logResults(leftIdx, rightIdx, diff, rprob);
		}
		// System.out.println("lIdx: " + leftIdx + ", rIdx: " + rightIdx + "\tmaxp: " + maxp.doubleValue() + "\tdiff: " + diff.doubleValue() + "\trprob: " + rprob.doubleValue());
		
		// System.out.println("pairedProbs:" + Util.nL() + Util.print2DArray(pairedProbs));
		
		return output;
	}

	public PPOutputHelix calculateCT(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs,
			double[][] distances, double alpha, double tau, int[] structure, boolean[][] canPair){
		
		CoFoldPosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new CoFoldPosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculate(insideProbs, outsideProbs, nucleotideProbs, distances, alpha, tau, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		PointRes[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		PointRes[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutputHelix output = null;
		
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutputHelix(new Helix(),null,PointRes.ZERO,PointRes.ZERO);
		} else {
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			PointRes[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			PointRes diff = diffs[leftIdx][rightIdx];
			//Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx, pairedProbs,diffs,canPair);
			output = new PPOutputHelix(helix, diffs, diff, rprob);
			
			logResults(leftIdx, rightIdx, diff, rprob);			
		}

		return output;
	}

	public PPOutputHelix calculateCTE(InsideOutsideProbabilities insideProbs,
			InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair){
		
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateE(insideProbs, outsideProbs, nucleotideProbs, structure, canPair);
		int leftIdx = posteriorProbabilities.getMaxLeftIdx();
		int rightIdx = posteriorProbabilities.getMaxRightIdx();
		PointRes[] unpairedProbs = posteriorProbabilities.getUnpairedProbs();
		PointRes[][] pairedProbs = posteriorProbabilities.getPairedProbs();
		PPOutputHelix output = null;
						
		if ((leftIdx<0)||(rightIdx<0)) {
			output = new PPOutputHelix(new Helix(),null,PointRes.ZERO,PointRes.ZERO);
		} else {
			boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, leftIdx, rightIdx);
			PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
			PointRes[][] diffs = posteriorProbabilitiesCalculator.getDiffs(pairedProbs, unpairedProbs, canPair);
			PointRes diff = diffs[leftIdx][rightIdx];
			//Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair);
			Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,pairedProbs,diffs,canPair);
			output = new PPOutputHelix(helix, diffs, diff, rprob);
			
			logResults(leftIdx, rightIdx, diff, rprob);
		}
		// System.out.println("lIdx: " + leftIdx + ", rIdx: " + rightIdx + "\tmaxp: " + maxp.doubleValue() + "\tdiff: " + diff.doubleValue() + "\trprob: " + rprob.doubleValue());
		
		// System.out.println("pairedProbs:" + Util.nL() + Util.print2DArray(pairedProbs));
		
		return output;
	}

	
	public PPOutputDouble calculateHybrid(InsideOutsideProbabilities insideProbs, InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
			double[][] distances, double weight, int[] structure, boolean[][] canPair, String alignmentPath, boolean[] blockedColumns) {
		
		int length = structure.length;
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		ShortDoublePosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateShort(insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
		double[] gUnpairedProbs = posteriorProbabilities.getUnpairedProbs();
		double[][] gPairedProbs = posteriorProbabilities.getPairedProbs();
		String structureString = dumpStructureWBC(structure,blockedColumns);
		double[][] vPairedProbs = 
				new Fasta2StockholmConvertingViennaRnaRunner(Settings.get().getViennaPath(), alignmentPath).getViennaProbabilities(structureString);
		// ProgramOutput.outMsg("vienna progs mtx:%n%s", Util.print2DArray(vPairedProbs));
		double[] vUnpairedProbs = getUnpairingProbs(vPairedProbs);
		double[] unpairedProbs = averageUnpairedProbs(gUnpairedProbs,vUnpairedProbs);
		double[][] pairedProbs = averagePairedProbs(gPairedProbs,vPairedProbs);
		double[][] vDiffs = getDiffs(vPairedProbs, vUnpairedProbs);
		double[][] diffs = getDiffs(pairedProbs, unpairedProbs);
		double[][] gDiffs = getDiffs(gPairedProbs, gUnpairedProbs);
		int leftIdx = -1; int rightIdx = -1;
		double maxp = 0;
		for (int j = 0; j<length; j++) {
			if (Constants.UnpairedBaseIdx == structure[j]) {
				for (int k = j+1; k<length; k++) {
					if (maxp<pairedProbs[j][k]) {
						leftIdx = j; rightIdx = k; 
						maxp = pairedProbs[j][k];
					}
				}
			}	
		}
		// double diff = diffs[leftIdx][rightIdx];
		double diff = Math.max(gDiffs[leftIdx][rightIdx], vDiffs[leftIdx][rightIdx]);
		double minDiff = Math.min(gDiffs[leftIdx][rightIdx], vDiffs[leftIdx][rightIdx]); 
		if ((minDiff<0)&&(diff<0.9)) {
			diff = -1;
		}
		
		Helix helix = new HelicesMaker().makeHelix(leftIdx, rightIdx, diffs, canPair);
		PPOutputDouble ppProbs = new PPOutputDouble(helix, diff);
		return ppProbs;
	}
	
	public PPOutputDouble calculateHybridE(InsideOutsideProbabilities insideProbs, InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
			int[] structure, boolean[][] canPair, String alignmentPath, boolean[] blockedColumns) {
		
		int length = structure.length;
		PosteriorProbabilitiesCalculator posteriorProbabilitiesCalculator = new PosteriorProbabilitiesCalculator(grammar);
		ShortDoublePosteriorProbabilities posteriorProbabilities = posteriorProbabilitiesCalculator.calculateShortE(insideProbs, outsideProbs, nucleotideProbs, structure, canPair);
		double[] gUnpairedProbs = posteriorProbabilities.getUnpairedProbs();
		double[][] gPairedProbs = posteriorProbabilities.getPairedProbs();
		String structureString = dumpStructureWBC(structure,blockedColumns);
		double[][] vPairedProbs = 
				//new OsCommandBasedViennaRnaRunner(SettingsFactory.getSettings().getViennaPath(), new SimpleFasta2StockholmConverter(alignmentPath)).getViennaProbabilities(alignmentPath, structureString);
				new Fasta2StockholmConvertingViennaRnaRunner(Settings.get().getViennaPath(), alignmentPath).getViennaProbabilities(structureString);
		double[] vUnpairedProbs = getUnpairingProbs(vPairedProbs);
		double[] unpairedProbs = averageUnpairedProbs(gUnpairedProbs,vUnpairedProbs);
		double[][] pairedProbs = averagePairedProbs(gPairedProbs,vPairedProbs);
		double[][] diffs = getDiffs(pairedProbs, unpairedProbs);
		int leftIdx = -1; int rightIdx = -1;
		double maxp = 0;
		for (int j = 0; j<length; j++) {
			if (Constants.UnpairedBaseIdx == structure[j]) {
				for (int k = j+1; k<length; k++) {
					if (maxp<pairedProbs[j][k]) {
						leftIdx = j; rightIdx = k; 
						maxp = pairedProbs[j][k];
					}
				}
			}	
		}
		double diff = diffs[leftIdx][rightIdx];
		Helix helix = new HelicesMaker().makeHelix(leftIdx, rightIdx, diffs, canPair);
		PPOutputDouble ppProbs = new PPOutputDouble(helix, diff);
		if(log.isDebugEnabled()){
			log.debug("gUnpairedProbs:%n{}", Util.print1DArray(gUnpairedProbs, 3, 5));
			log.debug("gPairedProbs:%n{}", Util.print2DArray(gPairedProbs, 5, 3));
			log.debug("vPairedProbs:%n{}", Util.print2DArray(vPairedProbs, 5, 3));
			log.debug("vUnpairedProbs:%n{}", Util.print1DArray(vUnpairedProbs, 3, 5));
			log.debug("unpairedProbs:%n{}", Util.print1DArray(unpairedProbs, 3, 5));
			log.debug("pairedProbs:%n{}", Util.print2DArray(pairedProbs, 5, 3));
			log.debug("diffs:%n{}", Util.print2DArray(diffs, 5, 3));
		}
		return ppProbs;
	}
	
	/**
	 * TODO: put the into right place 
	 * 
	 * @param structure
	 * @param blockedColumns
	 * @return
	 */
	private String dumpStructureWBC(int[] structure, boolean[] blockedColumns) {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<structure.length; i++) {
			char currentIdx = '+';
			if (structure[i] == -1) {
				if (blockedColumns[i]) {
					currentIdx = 'x';
				} else {
					currentIdx = '.';
				}
			} else if(structure[i] > i) {
				currentIdx = '(';
			} else if(structure[i] < i) {
				currentIdx = ')';
			} else {
				throw new InputMismatchException("Input was not an int array representing a structure");
			}	
			sb.append(currentIdx);
		}	
		String string = sb.toString();
		return string;
	}
	
	private double[][] averagePairedProbs(double[][] gPairedProbs, double[][] vPairedProbs) {
		int length = gPairedProbs.length;
		if (length != vPairedProbs.length) {
			throw new IllegalArgumentException("The two arrays must have the same dimensions");
		}
		double[][] pairedProbs = new double[length][length];
		for (int j = 0; j<length; j++) {
			for (int k = j+1; k<length; k++) {
				pairedProbs[j][k] = (gPairedProbs[j][k]+vPairedProbs[j][k])/2;
			}
		}
		return pairedProbs;
	}
	
	private double[] averageUnpairedProbs(double[] gUnpairedProbs, double[] vUnpairedProbs) {
		int length = gUnpairedProbs.length;
		if (length != vUnpairedProbs.length) {
			throw new IllegalArgumentException(String.format(
					"The two arrays must have the same dimensions. Got grammar len: %d, vienna len: %d",
					length, vUnpairedProbs.length ));
		}
		double[] unpairedProbs = new double[length];
		for (int j = 0; j<length; j++) {
			unpairedProbs[j] = (gUnpairedProbs[j]+vUnpairedProbs[j])/2;
		}
		return unpairedProbs;
	}
	
	private double[] getUnpairingProbs(double[][] pairingProbs) {
		double[] unpairingProbs = new double[pairingProbs.length];
		for (int j = 0; j<unpairingProbs.length; j++) {
			double tmp = 1.0;
			for (int k = 0; k<j; k++) {
				tmp = tmp - pairingProbs[k][j];
			}
			for (int k = j+1; k<pairingProbs.length; k++) {
				tmp = tmp - pairingProbs[j][k];
			}
			unpairingProbs[j] = tmp >= 0 ? tmp : 0;
//			if(tmp < 0){
//				throw new IllegalStateException(String.format(
//						"Got negative tmp: %g. j: %d. [%s]", tmp, j, Util.print1DArray(unpairingProbs, 3, 5)));
//			}
		}
		
		return unpairingProbs;
	}
	
	private double[][] getDiffs(double[][] pairingProbs, double[] unpairingProbs) {
		double[][] diffs = new double[pairingProbs.length][pairingProbs.length];
		for (int j = 0; j<pairingProbs.length; j++) {
			for (int k = j+1; k<pairingProbs.length; k++) {
				diffs[j][k] = pairingProbs[j][k]-(unpairingProbs[j]+unpairingProbs[k])/2;
			}
		}
		for(int i = 0; i < diffs.length; i++){
			for(int j = 0; j < diffs[i].length;j++){
				if(diffs[i][j]>1.001 || diffs[i][j]<-1.001){
					throw new IllegalStateException(String.format(
							"Got an expectation outside of allowed range [-1;+1] %.9g", diffs[i][j]));
				}
			}
		}
		return diffs;
	}
	
	private void logResults(int leftIdx, int rightIdx, PointRes diff, PointRes rprob) {
		if(log.isInfoEnabled()){
			log.info(String.format(Locale.UK, "lIdx: %3d; rIdx: %3d; diff: %7.4g; rprob: %7.4g",
					leftIdx, rightIdx, diff.doubleValue(), rprob.doubleValue()));
			//log.info(String.format(Locale.UK, "lIdx: %3d; rIdx: %3d; diff: %7.4g; rprob: %7.4g",
					//leftIdx, rightIdx, diff, rprob));
		}
	}


}
