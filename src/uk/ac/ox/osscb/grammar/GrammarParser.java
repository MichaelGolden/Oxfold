package uk.ac.ox.osscb.grammar;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Pattern;

import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.StochasticContextFreeGrammar;
import static java.util.regex.Pattern.CASE_INSENSITIVE;;

public class GrammarParser {
	
	public static HashMap<String, StochasticContextFreeGrammar> cached = new HashMap<String, StochasticContextFreeGrammar>();
	
	/**
	 * parses a grammar file and produces the corresponding grammar object to be used later
	 * a grammar file must contain the following lines:
	 * 1) !nonterminals A,B,C
	 * 			where A,B,C is the list of non-terminals (capital letters)
	 * 2) !terminals a,b,c
	 * 			where a,b,c is the list of terminals (lower case letters)
	 * 3) !rule A>BC|aDa>0.3|0.7
	 * 			where the non-terminal A can map to BC or aDa, with probabilities 0.3 and 0.7, respectively
	 * 			rules must be of one of the forms A->a, A->BC, A->aBb
	 * There must be non-terminals and terminals, and at least one rule. (If you're stupid and include more than one non-terminals line, only the last one is processed.)
	 * [In general, the idea is that lines starting with a ! should be lines containing information about the grammar, but the parser parses sloppy use of this notation.] 
	 */
	private static Pattern nonterminal = Pattern.compile("^!nonterminals ",CASE_INSENSITIVE);
	private static Pattern terminal = Pattern.compile("^!terminals ",CASE_INSENSITIVE);
	private static Pattern rule = Pattern.compile("^!rule ",CASE_INSENSITIVE);
	public static StochasticContextFreeGrammar parse(String gPath) {
		if(cached.containsKey(gPath))
		{
			return cached.get(gPath);
		}
		
		Scanner sc = null;
		try {
			sc = new Scanner(new File(gPath));
			List<ProductionRule> rules = new LinkedList<ProductionRule>();
			int TF,NTF,RF;
			TF = NTF = RF = 0; 
			char[] nts = new char[1],ts = new char[1];
			while (sc.hasNextLine()) {
				String line=sc.nextLine().replaceAll("\\n", "");
				if (nonterminal.matcher(line).find()) {
					NTF=1; line = line.replaceFirst("!nonterminals ", "");
					if (line.matches("^([A-Z]\\,)+[A-Z]$")) {
						String[] ntstmp = line.split("\\,");
						nts = new char[ntstmp.length];
						for (int j=0; j<ntstmp.length; j++) {
							nts[j]=ntstmp[j].charAt(0);
						}
					} else {
						throw new IllegalArgumentException("Wrong formatting of non-terminals string.");
					}
				} else if (terminal.matcher(line).find()) {
					TF=1; line = line.replaceFirst("!terminals ", "");
					if (line.matches("^([a-z]\\,)+[a-z]$")) {
						String[] tstmp = line.split("\\,");
						ts = new char[tstmp.length];
						for (int j=0; j<tstmp.length; j++) {
							ts[j]=tstmp[j].charAt(0);
						}
					} else {
						throw new IllegalArgumentException("Wrong formatting of terminals string.");
					}
				} else if (rule.matcher(line).find()) {
					RF = 1; line = line.replaceFirst("!rule ", "");
					if (line.matches("^[A-Z]\\>([a-zA-Z]+\\|)+[a-zA-Z]+\\>([01](\\.\\d+)?\\|)+[01](\\.\\d+)?$")) {
					// if (line.matches("^[A-Z]\\>([a-zA-Z]+\\|)+[a-zA-Z]+\\>([\\d\\.]+\\|)+[\\d\\.]+$")) {						
						String[] rtmp = line.split("\\>");
						String[] rtmpr = rtmp[1].split("\\|");
						String[] rtmpprob = rtmp[2].split("\\|");
						if (rtmpr.length==rtmpprob.length) {
							for (int j=0; j<rtmpr.length; j++) {
								char[] tmp = rtmp[0].toCharArray(); char tmp2 = tmp[0];
								char[] tmpr = rtmpr[j].toCharArray();
								String deleteme = rtmpprob[j];
								ProductionRule r = new ProductionRule(tmp2,tmpr,PointRes.valueOf(Double.valueOf(deleteme)));
								rules.add(r);
							}	
						} else {
							throw new IllegalArgumentException("Wrong formatting of rule string: number of possible right-hand sides does not match number of probabilities given.");
						}
					} else {
						throw new IllegalArgumentException("Wrong formatting of rule string.");
					}
				}	
			}
			if (0 == NTF) {
				throw new IllegalArgumentException("Idiot! You need to specify non-terminals.");
			} else if (0 == TF) {
				throw new IllegalArgumentException("Idiot! You need to specify terminals.");
			} else if (0 == RF) {
				throw new IllegalArgumentException("Idiot! You need to specify at least one rule.");
			}
			StochasticContextFreeGrammar g = new StochasticContextFreeGrammar(nts,ts,rules);
			cached.put(gPath, g);
			return g;
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",gPath));
		}finally{
			if(null != sc){
				sc.close();
			}
		}
	}
}
