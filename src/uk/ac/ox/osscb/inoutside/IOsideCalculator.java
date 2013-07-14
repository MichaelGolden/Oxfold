package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PossiblePairFinder;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;

public interface IOsideCalculator {

	/**
	 * @param pairingProbs
	 * @param distances
	 * @param weight
	 * @param structure - the pairing structure to be calculated
	 * @param canPair output of {@link PossiblePairFinder#canPair(int[])
	 * @return
	 */
	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair);
	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair);

	public InsideOutsideProbabilities insideE(
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair);
	
	/**
	 * 
	 * @param insideProbs - output of the {@link InsideOutsideCalculator#inside(PairingProbabilities, int[][], double)}
	 * function
	 * @param pairingProbs - 'n' in pseudocode
	 * @param distances - 'd'
	 * @param structure
	 * @param weight
	 * @return
	 */
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[][] distances,
			double weight, int[] structure, boolean[][] canPair);
	
	public InsideOutsideProbabilities outside(InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double weight, int[] structure, boolean[][] canPair);

	public InsideOutsideProbabilities outsideE(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair);
	

}