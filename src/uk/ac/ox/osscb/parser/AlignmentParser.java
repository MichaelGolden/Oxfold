package uk.ac.ox.osscb.parser;

import uk.ac.ox.osscb.Alphabet;

public interface AlignmentParser {
	
	/**
	 * @return alignments array, nucleotides plus gaps in form of dots 
	 * or other symbol (TODO: to consider).
	 */
	String[] parse(String alignmentFile);
	
	String[] parseEvolutionary(String alignmentFile, Alphabet alphabet);
}
