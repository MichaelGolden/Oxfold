package uk.ac.ox.osscb;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;

public interface DefaultParametersProvider {
	
	/**
	 * 
	 * @return 16 expected. Entered by user ideally.
	 */
	NucleotideProbsPrecise getDefaultBasingProbabilities();

}
