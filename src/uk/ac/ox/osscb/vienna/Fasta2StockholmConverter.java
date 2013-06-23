package uk.ac.ox.osscb.vienna;

/**
 * Used to convert fasta to a stockholm formats.
 * 
 * @author Vladimir
 *
 */
public interface Fasta2StockholmConverter {
	
	/**
	 * 
	 * @param outputFile if null, the resulting string is convertion
	 * in a text form only.
	 * @return text representation of the convertion.
	 */
	String convert(String outputFile);

}
