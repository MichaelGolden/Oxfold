package uk.ac.ox.osscb;

import uk.ac.ox.osscb.domain.NucleotideProbsDouble;

public class AlignmentConverter {

	/**
	 * converts alignment to array of integers for easy processing
	 * (the converter is quite lenient: any non-alphabet characters are parsed as gaps)
	 * (the gap character is assumed to be the first character of the alphabet)
	 * TODO do something more sensible about the above	
	 * @param alin - alignment of sequences; input
	 * @param param - parameter from parameter parser
	 * @return alignment in format for other stuff
	 */
	public int[][] convert(String[] alin, NucleotideProbsDouble param) {
		int dim = alin[0].length();
		int seqnum = alin.length;
		String alph = param.getAlphabetStr();
		int[][] alignment = new int[seqnum][dim];
		for (int j=0; j<seqnum; j++) {
			if (alin[j].length() != dim) {
				throw new IllegalArgumentException("You should really try to provide proper alignments: the sequences " +
						"in your alignment do not all have the same length.");
			}
			for (int k=0; k<dim; k++) {
				int tmp = alph.indexOf(alin[j].charAt(k)); 
				if (tmp>0) {
					alignment[j][k] = tmp;
				} else {
					alignment[j][k] = 0;
				}
			}
		}
		return alignment;
	}
}
