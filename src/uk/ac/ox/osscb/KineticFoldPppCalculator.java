package uk.ac.ox.osscb;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputHelixInternalResult;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResultDouble;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.ShortDoublePosteriorProbabilities;

/**
 * Base interface for high-level calculation of
 * {@link PPOutput}.
 * Implementations use either weight-aware or weight-less approaches.
 * 
 * @author Vladimir
 */
public interface KineticFoldPppCalculator {


	PPOutput calculateDynamicPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair);
	
	PPOutputHelixInternalResult calculateCotranscriptionalPpOutput(NucleotideProbsPrecise alignmentProbs, 
			int[] structure, boolean[][] canPair, PosteriorProbabilities postProbs);
	
	PPOutputInternalResultDouble calculateHybridPpOutput (NucleotideProbsPrecise alignmentProbs, int[] structure, 
			String alignmentPath, ShortDoublePosteriorProbabilities postProbs, boolean[] blockedColumns);

	PPOutputInternalResult2 calculatePpOutputInternalResult2(
			NucleotideProbsPrecise alignmentProbs, int[] structure,
			PosteriorProbabilities postProbs);

	PPOutput calculatePpOutput(NucleotideProbsPrecise alignmentProbs,
			int[] structure);


}
