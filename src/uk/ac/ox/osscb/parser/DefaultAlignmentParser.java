package uk.ac.ox.osscb.parser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.regex.Pattern;

import uk.ac.ox.osscb.Alphabet;

public class DefaultAlignmentParser implements AlignmentParser {

	private static Pattern newseq = Pattern.compile("^\\>");

	public String[] parse(String alignmentFile) {
		Scanner sc = null;
		try {
			sc = new Scanner(new File(alignmentFile));
			int seqnum = 0;
			String line = "";
			while (sc.hasNextLine()) {
				line = sc.nextLine();
				if (newseq.matcher(line).find()) {
					seqnum++;
				}
			}
			if ((0==seqnum)||(0==line.indexOf('>'))) {
				throw new IllegalArgumentException("No sequences in file (FASTA Format?) or last sequence missing.");
			}
			sc.close();
			Scanner scn = new Scanner(new File(alignmentFile));
			String[] align = new String[seqnum];
			boolean nextHasSeq = false;
			int found = 0;
			while (scn.hasNextLine()) {
				line = scn.nextLine();
				if (newseq.matcher(line).find()) {
					if (nextHasSeq) {
						throw new IllegalArgumentException("Missing sequence? Unable to parse.");
					} else {
						nextHasSeq = true;
					}
				} else if (nextHasSeq) {
					line = line.replaceAll("\\n", "");
					align[found] = line;
					found++;
					nextHasSeq = false;
				}
			}
			scn.close();
			return align;
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",alignmentFile));
		} finally {
			if (null != sc) {
				sc.close();
			}
		}
	}
	
	public String[] parseEvolutionary(String alignmentFile, Alphabet alphabet) {
		Scanner sc = null;
		try {
			sc = new Scanner(new File(alignmentFile));
			int seqnum = 0;
			String line = "";
			while (sc.hasNextLine()) {
				line = sc.nextLine();
				if (newseq.matcher(line).find()) {
					seqnum++;
				}
			}
			if ((0==seqnum)||(0==line.indexOf('>'))) {
				throw new IllegalArgumentException("No sequences in file (FASTA Format?) or last sequence missing.");
			}
			sc.close();
			Scanner scn = new Scanner(new File(alignmentFile));
			String[] align = new String[seqnum];
			boolean nextHasSeq = false;
			int found = 0;
			while (scn.hasNextLine()) {
				line = scn.nextLine();
				if (newseq.matcher(line).find()) {
					if (nextHasSeq) {
						throw new IllegalArgumentException("Missing sequence? Unable to parse.");
					} else {
						nextHasSeq = true;
					}
				} else if (nextHasSeq) {
					line = line.replaceAll("\\n", "");
					String seq = "";
					String[] tmp = line.split("");
					for (int j = 1; j<tmp.length; j++) {
						if (!alphabet.getSynonyms().containsKey(tmp[j])) {
							tmp[j] = alphabet.getNames()[0];
						}
						seq = seq.concat(tmp[j]);
					}
					align[found] = seq;
					found++;
					nextHasSeq = false;
				}
			}
			scn.close();
			return align;
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",alignmentFile));
		} finally {
			if (null != sc) {
				sc.close();
			}
		}
	}
	
	public static int calculateAlignmentLength(String alingmentFile){
		
		Scanner sc = null;
		try {
			
			for(sc = new Scanner(new File(alingmentFile));
					sc.hasNext();){

				String nextLine = sc.next().trim();
				if(0 == nextLine.length() || newseq.matcher(nextLine).matches()){
					// we've hit 'new sequence header' or an empty line, skip it in any case
					continue;
				}
				return nextLine.length();
			}
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("cannot find file named: %s", alingmentFile));
		} finally {
			if(null != sc)
				sc.close();
		}
		
		return 0;
	}
}
