package uk.ac.ox.osscb.vienna;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

import uk.ac.ox.osscb.Util;
import uk.ac.ox.osscb.util.TextUtil;

/**
 * Assumes the file is in the following format:
 * 
 * .......<<<<<.........>>>>>.............
 * >CP000886.1/4403693-4403730
 * AAGCGUAUUGGUAGCAG.UAAGCCAAGGGCGGUAGCGUU
 * >CP000826.1/3180995-3180958
 * UAGCGCAUUGGGAGCUG.UAACUCAAGGGCGGUAGCGUG
 * >AE014075.1/2424993-2424955
 * CAGUGUAUUGGUAGCUAAAAAGCCAGGGGCGGUAGCGUG
 * 
 * The 1st line is optional and can be omitted.
 * It is omitted during conversion.
 * Converter is rather strict.
 * 
 * @author Vladimir
 *
 */
public class SimpleFasta2StockholmConverter implements Fasta2StockholmConverter {

	private static final Pattern zeroLine = Pattern.compile("^[.<>]+$");
	private static final Pattern firstLine = Pattern.compile(">[A-Z][.A-Z0-9/-]+.+$");
	private static final Pattern secondLine = Pattern.compile("^[.ACDGKMNRSUVWY]+$");
	
	private String fastaFile;
	
	public SimpleFasta2StockholmConverter(String fastaFile) {
		super();
		this.fastaFile = fastaFile;
	}

	/**
	 * 
	 * @param outputFile if null, the resulting string is convertion
	 * in a text form only.
	 * 
	 * @return
	 */
	@Override
	public String convert(String outputFile) {
		Util.assertCanReadFile(this.fastaFile);
		
		String converted;
		Scanner sc = null;
		try {
			sc = new Scanner(new File(this.fastaFile));
			converted = readAndConvert(sc);
		} catch (FileNotFoundException ex) {
			// in theory we should never get here
			throw new IllegalArgumentException(String.format("file: '%s' cannot be found", this.fastaFile));
		}finally{
			if(null != sc){
				sc.close();
				sc = null;
			}
		}
		
		if(null != outputFile){
			TextUtil.writeTextFile(outputFile, converted);
		}
		
		return converted;
	}

	private String readAndConvert(Scanner sc){

		String seqName = null;
		String seqItself = null;

		StringBuilder sb = new StringBuilder();
		// we expect 1st line to be now. 2nd line - otherwise
		// we could expect either zero line (.<<..>> ...)
		boolean isFirstLine = true;
		for(int lineNo = 0;sc.hasNext(); lineNo++){
			
			String line = sc.next();
			if(0 == lineNo){
				if(this.zeroLine.matcher(line).matches()){
					// ok, skip the line and continue iterating
					continue;
				}
			}
			
			if(isFirstLine){
				// we want to have the 1st line now
				if(firstLine.matcher(line).matches()){
					seqName = line;
					isFirstLine = false;
					continue;
				} /*else if(line.trim().equals("//")){
					// we're done. check for correct state:
					if(!StringUtils.isBlank(seqName) || !StringUtils.isBlank(seqItself))
						throw new IllegalStateException(String.format(
								"wrong state: EOF while either name: '%s' or seq: '%s' has text", seqName, seqItself));
					// check there are no more lines left:
					if(sc.hasNext()){
						
					}
					break;
				}*/else{
					throw new IllegalStateException(String.format(
							"Cannot parse file: %s. Want first line, got non-matched: '%s'. Stockholm so far:%n%s",
							this.fastaFile, line, sb));
				}
			}else{
				
				if(secondLine.matcher(line).matches()){
					seqItself = line;
					isFirstLine = true;
				}else{
					throw new IllegalStateException(String.format(
							"Cannot parse file: %s. Want first line, got non-matched: '%s'. Stockholm so far:%n%s",
							this.fastaFile, line, sb));
				}
			}
			assertBothHaveText(seqName, seqItself);
			
			sb.append(Util.nL()).append(seqName).append(" ").append(seqItself);
			// reset seq data collected so far
			seqName = null;
			seqItself = null;
		}
		
		assertBothAreBlank(seqName, seqItself);
		sb.insert(0, "# STOCKHOLM 1.0");
		sb.append(Util.nL()).append("//");
		return sb.toString();
	}

	private void assertBothHaveText(String seqName, String seqItself) {
		if(StringUtils.isBlank(seqName))
			throw new IllegalStateException("attempt to add sequence with an empty name");
		if(StringUtils.isBlank(seqItself))
			throw new IllegalStateException("attempt to add an empty sequence");
	}
	
	private void assertBothAreBlank(String seqName, String seqItself) {
		if(!StringUtils.isBlank(seqName) || !StringUtils.isBlank(seqItself))
			throw new IllegalStateException(String.format(
					"wrong state. Expected name and sequence itself be blank. In fact it was: '%s' and '%s' correspondingly", 
					seqName, seqItself));
	}
}
