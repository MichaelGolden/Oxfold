package uk.ac.ox.osscb.analysis;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.HashMap;

/**
 *
 * @author Michael Golden <michaelgolden0@gmail.com>
 */
public class AmbiguityCodes {

    static HashMap<String, String> codes = new HashMap<String, String>();

    public AmbiguityCodes() {
        codes.put("A", "A");
        codes.put("C", "C");
        codes.put("G", "G");
        codes.put("T", "T");
        codes.put("U", "T"); 
        codes.put("Y", "CT");
        codes.put("R", "AG");
        codes.put("W", "AT");
        codes.put("S", "GC");
        codes.put("K", "TG");
        codes.put("M", "CA");
        codes.put("D", "AGT");
        codes.put("V", "ACG");
        codes.put("H", "ACT");
        codes.put("B", "CGT");
        codes.put("N", "ACGT");
        codes.put("X", "ACGT");
        codes.put("-", "-");
    }

    public double[] getBaseScores(String ambiguityCode) {
        double[] scores = new double[5];
        String standardBases = codes.get(ambiguityCode);
        if (standardBases != null) {
            for (int i = 0; i < standardBases.length(); i++) {
                switch (standardBases.charAt(i)) {
                    case 'A':
                        scores[0] += 1;
                        break;
                    case 'C':
                        scores[1] += 1;
                        break;
                    case 'G':
                        scores[2] += 1;
                        break;
                    case 'T':
                        scores[3] += 1;
                        break;
                    case '-':
                        scores[4] += 1;
                        break;
                    default:
                        break;
                }
            }
            for (int i = 0; i < scores.length; i++) {
                scores[i] /= ((double) standardBases.length());
            }
        } else {
            System.out.println(ambiguityCode);
        }

        return scores;
    }
    
    public char[] getAmbigChars(char ambiguityCode){
    	String standardBases = codes.get(Character.toString(ambiguityCode));
    	return standardBases.toCharArray(); 
    }
}
