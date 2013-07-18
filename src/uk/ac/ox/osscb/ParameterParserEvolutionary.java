package uk.ac.ox.osscb;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Pattern;

public class ParameterParserEvolutionary {
	public static HashMap<String, EvolutionaryParameters> cached = new HashMap<String, EvolutionaryParameters>();
	
	
	/**
	 * parses parameter file
	 */
	private static Pattern alph = Pattern.compile("^!alphabet ");
	private static Pattern unpaired = Pattern.compile("^!unpaired");
	private static Pattern paired = Pattern.compile("^!paired");
	private static Pattern special = Pattern.compile("^!special ");
	private static Pattern testend = Pattern.compile("^\\?");
	
	public static EvolutionaryParameters parse(String pPath) {
		if(cached.containsKey(pPath))
		{
			return cached.get(pPath);
		}
		
		int AF = 0,UPPF = 0,PPF =0,SF=0;
		Letter[] letters = new Letter[1];
		Alphabet alphabet = null; Alphabet pAlphabet = null;
		double[][] unpairedR = null; double[][] unpairedD = null; double[][] unpairedRI = null; double[] unpairedprior = null;
		double[][] pairedR = null; double[][] pairedD = null; double[][] pairedRI = null; double[] pairedprior = null;
		String alphstring = "";
		try {
			Scanner sc = new Scanner(new File(pPath));
			while (sc.hasNextLine()) {
				String line = sc.nextLine().replaceAll("\n", "");
				if (alph.matcher(line).find()) {
					AF = 1; line = line.replaceFirst("!alphabet ","");
					if (line.matches("^(\\w\\,)*[\\w]$")) {
						String[] tmp = line.split("\\,");
						letters = new Letter[tmp.length];
						for (int j = 0; j<tmp.length; j++) {
							letters[j] = new Letter(tmp[j],new String[]{tmp[j]});
							alphstring = alphstring.concat(tmp[j]);
						}
					} else {
						throw new IllegalArgumentException("Alphabet string has wrong format.");
					}
				} else if (special.matcher(line).find()) {
					SF = 1;
					if (AF == 0) {
						throw new IllegalArgumentException("Need to specify alphabet first.");						
					} else {
						SF = 1; line = line.replaceFirst("!special ","");
						String[] specialRules = line.split("\\;");
						for (int j = 0; j<specialRules.length; j++) {
							String[] tmp = specialRules[j].split(" ");
							int idx = alphstring.indexOf(tmp[0]);
							letters[idx] = new Letter(tmp[0],Arrays.copyOfRange(tmp, 1, tmp.length));
						}
						alphabet = new Alphabet(letters);
						int length = letters.length;
						Letter[] pLetters = new Letter[length*length];
						for (int j = 0; j<length; j++) {
							for (int k = 0; k<length; k++) {
								String tmpname = letters[j].getName().concat(letters[k].getName());
								String[] synj = letters[j].getSynonyms();
								String[] synk = letters[k].getSynonyms();
								String[] tmpsyn = new String[synj.length*synk.length];
								for (int l = 0; l<synj.length; l++) {
									for (int m = 0; m<synk.length; m++) {
										tmpsyn[l*synk.length + m] = synj[l].concat(synk[m]);  
									}
								}
								pLetters[j*length+k] = new Letter(tmpname,tmpsyn);
							}
						}
						pAlphabet = new Alphabet(pLetters);
					}		
				} else if (unpaired.matcher(line).find()) {
					if (AF*SF == 0) {
						throw new IllegalArgumentException("Need to specify alphabet and special characters before specifying unpaired probabilities.");
					} else {
						UPPF = 1; line = sc.nextLine().replaceAll("\n", "");
						String[] tmp = line.split(" ");
						int length = alphabet.getStandardNames().length;
						unpairedprior = new double[length];
						unpairedR = new double[length][length];
						unpairedD = new double[length][length];
						unpairedRI = new double[length][length];
						if (tmp.length != length) {
							throw new IllegalArgumentException("Number of probabilities specified does not match length of alphabet.");
						}
						for (int j = 0; j<tmp.length; j++) {
							unpairedprior[j] = Double.valueOf(tmp[j]);
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								unpairedR[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								unpairedD[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								unpairedRI[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine().replaceAll("\n", "");
						if (!testend.matcher(line).find()) {
							throw new IllegalArgumentException("Unpaired probabilities could not be parsed.");
						}
					}
				} else if (paired.matcher(line).find()) {
					if (AF*SF == 0) {
						throw new IllegalArgumentException("Need to specify alphabet and special characters before specifying paired probabilities.");
					} else {
						PPF = 1; line = sc.nextLine().replaceAll("\n", "");
						String[] tmp = line.split(" ");
						int length = pAlphabet.getStandardNames().length;
						pairedprior = new double[length];
						pairedR = new double[length][length];
						pairedD = new double[length][length];
						pairedRI = new double[length][length];
						if (tmp.length != length) {
							throw new IllegalArgumentException("Number of probabilities specified does not match length of alphabet.");
						}
						for (int j = 0; j<tmp.length; j++) {
							pairedprior[j] = Double.valueOf(tmp[j]);
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								pairedR[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								pairedD[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine();
						for (int j = 0; j<length; j++) {
							line = sc.nextLine().replaceAll("\n", "");
							tmp = line.split(" ");
							for (int k = 0; k<length; k++) {
								pairedRI[j][k] = Double.valueOf(tmp[k]);
							}
						}
						line = sc.nextLine().replaceAll("\n", "");
						if (!testend.matcher(line).find()) {
							throw new IllegalArgumentException("Paired probabilities could not be parsed.");
						}
					}
				}
			}
			if (AF*UPPF*PPF*SF == 0) {
				throw new IllegalArgumentException("Need to specify alphabet, special characters, unpaired probabilities and paired probabilities");
			} 
			EvolutionaryParameters param = new EvolutionaryParameters(alphabet, new Qmatrix(unpairedR,unpairedD,unpairedRI,unpairedprior), pAlphabet, new Qmatrix(pairedR,pairedD,pairedRI,pairedprior));
			cached.put(pPath, param);
			return param;
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",pPath));
		}	
	}
}
