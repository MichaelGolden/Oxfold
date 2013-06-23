package uk.ac.ox.osscb.vienna;

public interface ViennaRunner {

	/**
	 * @param structureConstraint - the string to be passed to RNAalifold to stdin.
	 * See RNAalifold <i>-C, --constraint</i> option description.
	 * @return
	 */
	double[][] getViennaProbabilities(String constrainingStructure);
	
}
