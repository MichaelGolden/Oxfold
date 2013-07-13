package uk.ac.ox.osscb.analysis;
import java.util.Arrays;

/**
 * Based on Moulton et al.
 *
 * @author Michael Golden <michaelgolden0@gmail.com>
 */
public class MountainMetrics {

    public static void main(String[] args) {
        /*
         * int[] s3 =
         * RNAFoldingTools.getPairedSitesFromDotBracketString("(......)"); int[]
         * s4 = RNAFoldingTools.getPairedSitesFromDotBracketString("..(..)..");
         * System.out.println(MountainMetrics.calculateMountainDistance(s3,
         * MountainMetrics.getStructureZero(8), 1));
         * System.out.println(MountainMetrics.calculateMountainDistance(s4,
         * MountainMetrics.getStructureZero(8), 1));
         */

        int[] s1 = RNAFoldingTools.getPairedSitesFromDotBracketString("...(..).......");
        int[] sz = RNAFoldingTools.getPairedSitesFromDotBracketString("..............");
        System.out.println(calculateWeightedMountainDistance(s1, sz));
    }
    
    public static double calculateWeightedMountainSimilarity(int[] pairedSites1, int[] pairedSites2)
    {
    	return 1 - calculateNormalizedWeightedMountainDistance(pairedSites1, pairedSites2);
    }

    public static double calculateNormalizedWeightedMountainDistance(int[] pairedSites1, int[] pairedSites2) {
        return calculateWeightedMountainDistance(pairedSites1, pairedSites2) / calculateWeightedMountainDiameter(pairedSites1.length);
    }

    public static double calculateWeightedMountainDistance(int[] pairedSites1, int[] pairedSites2) {
        double[] f1 = getMountainVector(pairedSites1, true);
        double[] f2 = getMountainVector(pairedSites2, true);
        double d = 0;
        for (int i = 0; i < f1.length; i++) {
            d += Math.abs(f1[i] - f2[i]);
        }
        return d;
    }

    public static double calculateWeightedMountainDiameter(int length) {
        return calculateWeightedMountainDistance(getStructureStar(length), getStructureZero(length));
    }

    public static double calculateMountainDistance(int[] pairedSites1, int[] pairedSites2, int p) {
        double[] f1 = getMountainVector(pairedSites1, false);
        double[] f2 = getMountainVector(pairedSites2, false);
        double d = 0;
        for (int i = 0; i < f1.length; i++) {
            d += Math.pow(Math.abs(f1[i] - f2[i]), p);
        }
        return Math.pow(d, 1 / p);
    }

    /**
     * Optimised mountain vector calculation.
     *
     * @param pairedSites
     * @param weighted
     * @return
     */
    public static double[] getMountainVector(int[] pairedSites, boolean weighted) {
        double[] f1 = new double[pairedSites.length];
        if (weighted) {
            if (pairedSites[0] != 0) {
                f1[0] += 1 / ((double) (pairedSites[0] - 1 - 0));
            }

            for (int i = 1; i < pairedSites.length; i++) {
                f1[i] = f1[i - 1];
                if (pairedSites[i] != 0) {
                    f1[i] += 1 / ((double) (pairedSites[i] - 1 - i));
                }
            }
            return f1;
        } else {
            if (pairedSites[0] != 0) {
                if (pairedSites[0] > 0) {
                    f1[0] += 1;
                } else {
                    f1[0] -= 1;
                }
            }

            for (int i = 1; i < pairedSites.length; i++) {
                f1[i] = f1[i - 1];
                if (pairedSites[i] != 0) {
                    if (pairedSites[i] > i) {
                        f1[i] += 1;
                    } else {
                        f1[i] -= 1;
                    }
                }
            }
            return f1;
        }
    }

    public static double[] getMountainVector2(int[] pairedSites, boolean weighted) {
        String s1 = RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites);
        double[] f1 = new double[pairedSites.length];

        double w = 1;
        if (weighted) {
            w = 1 / (double) Math.abs(pairedSites[0] - 1 - 0);
        }
        if (s1.charAt(0) == '(') {
            f1[0] += w;
        } else if (s1.charAt(0) == ')') {
            f1[0] -= w;
        }

        for (int i = 1; i < pairedSites.length; i++) {
            if (weighted) {
                w = 1 / (double) Math.abs(pairedSites[i] - 1 - i);
            }

            f1[i] = f1[i - 1];
            if (s1.charAt(i) == '(') {
                f1[i] += w;
            } else if (s1.charAt(i) == ')') {
                f1[i] -= w;
            }
        }
        return f1;
    }

    public static double[] compareMountainSlopes(int[] pairedSites1, int[] pairedSites2, boolean weighted) {
        double[] mountain1 = MountainMetrics.getMountainVector(pairedSites1, weighted);
        double[] mountain2 = MountainMetrics.getMountainVector(pairedSites2, weighted);

        int windowSize = 10;
        double[] slope1 = getMountainSlopeVector(pairedSites1, weighted, windowSize);
        double[] slope2 = getMountainSlopeVector(pairedSites2, weighted, windowSize);

        for (int i = 0; i < slope1.length; i++) {
            double angle1 = Math.acos(slope1[i]);
            double angle2 = Math.acos(slope2[i]);
            double sim = 1 - (Math.abs(angle1 - angle2) / Math.PI);


            System.out.println(i+"\t"+mountain1[i]+"\t"+mountain2[i]+"\t"+slope1[i]+"\t"+slope2[i]+"\t"+angle1 + "\t" + angle2 + "\t" + sim);
        }

        return null;
    }

    public static double[] getMountainSlopeVector(int[] pairedSites, boolean weighted, int windowSize) {
        double[] mountain = MountainMetrics.getMountainVector(pairedSites, weighted);
        double[] slope = new double[mountain.length];
        for (int i = windowSize ; i < mountain.length; i++) {
            slope[i] = (mountain[i] - mountain[i - windowSize])/((double)windowSize);
        }
        return slope;
    }

    public static int[] getStructureStar(int length) {
        String dotBracket;
        if (length % 2 == 0) // even
        {
            dotBracket = nChars('(', (length - 2) / 2) + ".." + nChars(')', (length - 2) / 2);
        } else // odd
        {
            dotBracket = nChars('(', (length - 1) / 2) + "." + nChars(')', (length - 1) / 2);
        }
        return RNAFoldingTools.getPairedSitesFromDotBracketString(dotBracket);
    }

    public static int[] getStructureZero(int length) {
        int[] structureZero = new int[length];
        Arrays.fill(structureZero, 0); // probably not necessary
        return structureZero;
    }

    public static String nChars(char c, int length) {
        String ret = "";
        for (int i = 0; i < length; i++) {
            ret += c;
        }
        return ret;
    }
}
