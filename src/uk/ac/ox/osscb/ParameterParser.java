package uk.ac.ox.osscb;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.regex.Pattern;

import uk.ac.ox.osscb.domain.NucleotideProbsDouble;

public class ParameterParser {
	/**
	 * parses parameter file
	 */
	private static Pattern alph = Pattern.compile("^!alphabet ");
	private static Pattern unpaired = Pattern.compile("^!unpaired ");
	private static Pattern paired = Pattern.compile("^!paired ");
	private static Pattern special = Pattern.compile("^!special ");
	public NucleotideProbsDouble parse(String pPath) {
		int AF = 0,UPPF = 0,PPF =0,SF=0;
		char[] alphabet = new char[1];
		double[][] PairedProbs = new double[1][1];
		double[] UnpairedProbs = new double[1];
		double[][] PairedProbsSpecial = new double[1][1];
		double[] UnpairedProbsSpecial = new double[1];
		String alphstring = "";
		try {
			Scanner sc = new Scanner(new File(pPath));
			while (sc.hasNextLine()) {
				String line = sc.nextLine().replaceAll("\n", "");
				if (alph.matcher(line).find()) {
					AF = 1; line = line.replaceFirst("!alphabet ","");
					if (line.matches("^([^,]\\,)*[^,]$")) {
						String[] tmp = line.split("\\,");
						alphabet = new char[tmp.length];
						for (int j = 0; j<tmp.length; j++) {
							alphabet[j] = tmp[j].charAt(0);
						}
						alphstring = new String(alphabet);
					} else {
						throw new IllegalArgumentException("Alphabet string has wrong format.");
					}
				} else if (unpaired.matcher(line).find()) {
					UPPF = 1; line = line.replaceFirst("!unpaired ", "");
					if (line.matches("^([01](\\.\\d+)?\\,)*[01](\\.\\d+)?$")) {
						String[] tmp = line.split("\\,");
						UnpairedProbs = new double[tmp.length];
						for (int j = 0; j<tmp.length; j++) {
							UnpairedProbs[j] = Double.valueOf(tmp[j]);
						}
					} else  {
						throw new IllegalArgumentException("Wrong formatting of Unpaired Probabilities string.");
					}
				} else if (paired.matcher(line).find()) {
					PPF = 1; line = line.replaceFirst("!paired ","");
					if (line.matches("^(([01](\\.\\d+)?\\,)*[01](\\.\\d+)?\\;)*([01](\\.\\d+)?\\,)*[01](\\.\\d+)?$")) {
						String[] tmp = line.split("\\;");
						String[] ttmp = tmp[0].split("\\,");
						int l = ttmp.length;
						PairedProbs = new double[tmp.length][l];
						for (int k=0;k<l;k++) {
							PairedProbs[0][k] = Double.valueOf(ttmp[k]);
						}
						for (int j=1;j<tmp.length; j++) {
							ttmp = tmp[j].split("\\,");
							if (l == ttmp.length) {
								for (int k=0;k<l;k++) {
									PairedProbs[j][k] = Double.valueOf(ttmp[k]);
								}
							} else {
								throw new IllegalArgumentException("Badly formed Paired Probabilities string: missed some probabilities?");
							}
						}
					} else {
						throw new IllegalArgumentException("Wrong formatting of Paired Probabilities string.");
					}
				} else if (special.matcher(line).find()) {
					if (AF*UPPF*PPF == 0) {
						throw new IllegalArgumentException("Need to specify everything else first.");						
					} else {
						SF=1; line = line.replaceFirst("!special ", "");
						String[] tmp = line.split("\\;");
						PairedProbsSpecial = new double[alphabet.length][alphabet.length];
						UnpairedProbsSpecial = new double[alphabet.length];
						for (int j = 0; j<UnpairedProbs.length; j++) {
							UnpairedProbsSpecial[j] = UnpairedProbs[j];
							for (int k = 0; k<UnpairedProbs.length; k++) {
								PairedProbsSpecial[j][k]=PairedProbs[j][k];
							}
						}
						for (int j = 0; j<tmp.length; j++) {
							String[] tmpp = tmp[j].split(" ");
							double sumup = 0; 							
							for (int k = 1; k<tmpp.length; k++) {
								sumup = sumup + UnpairedProbs[alphstring.indexOf(tmpp[k])];
							}
							UnpairedProbsSpecial[alphstring.indexOf(tmpp[0])] = sumup/(tmpp.length - 1);
							for (int k = 0; k<UnpairedProbs.length; k++) {
								double sumpl = 0,sumpr = 0;
								for (int l = 1; l<tmpp.length; l++) {
									sumpl = sumpl + PairedProbs[alphstring.indexOf(tmpp[l])][k];
									sumpr = sumpr + PairedProbs[k][alphstring.indexOf(tmpp[l])];
								}
								PairedProbsSpecial[k][alphstring.indexOf(tmpp[0])] = sumpr/(tmpp.length-1);
								PairedProbsSpecial[alphstring.indexOf(tmpp[0])][k] = sumpl/(tmpp.length-1);
							}
						}
						for (int j = 0; j<tmp.length; j++) {
							String[] tmpp = tmp[j].split(" ");
							for (int k = PairedProbs.length; k<alphabet.length; k++) {
								double sumpl = 0, sumpr = 0;
								for (int l = 1; l<tmpp.length; l++) {
									sumpl = sumpl + PairedProbsSpecial[alphstring.indexOf(tmpp[l])][k];
									sumpr = sumpr + PairedProbsSpecial[k][alphstring.indexOf(tmpp[l])];
								}
								PairedProbsSpecial[alphstring.indexOf(tmpp[0])][k] = sumpl/(tmpp.length-1);
								PairedProbsSpecial[k][alphstring.indexOf(tmpp[0])] = sumpr/(tmpp.length-1);
							}
						}	
					}
				}
			}
			if (SF == 0) {
				PairedProbsSpecial = PairedProbs; UnpairedProbsSpecial=UnpairedProbs;
			}
			if (AF*UPPF*PPF == 0) {
				throw new IllegalArgumentException("Need to specify alphabet, unpaired probabilities and paired probabilities.");
			} else if ((alphabet.length != UnpairedProbsSpecial.length)||(alphabet.length != PairedProbsSpecial.length)) {
				throw new IllegalArgumentException("Length of alphabet must match dimensions of unpaired and paired probabilities.");
			}
			NucleotideProbsDouble param = new NucleotideProbsDouble(PairedProbsSpecial,UnpairedProbsSpecial,alphabet);
			return param;
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",pPath));
		}	
	}
}
