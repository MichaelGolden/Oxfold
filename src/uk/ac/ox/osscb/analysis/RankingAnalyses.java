package uk.ac.ox.osscb.analysis;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Michael
 */
public class RankingAnalyses {

    public static double[] toDoubleArray(ArrayList<Double> array) {
        double[] values = new double[array.size()];
        for (int i = 0; i < array.size(); i++) {
            values[i] = array.get(i);
        }
        return values;
    }
    
    public static double [] getArray(List<Double> values)
    {
        double [] ret = new double[values.size()];
        for(int i = 0 ; i < ret.length ; i++)
        {
            ret[i] = values.get(i);
        }
        return ret;
    }

    public static double getAverage(ArrayList<Double> list) {
        double sum = 0;
        for (int i = 0; i < list.size(); i++) {
            sum += list.get(i);
        }

        return sum / (double) list.size();
    }

    public static double getAverage(double[] list) {
        double sum = 0;
        for (int i = 0; i < list.length; i++) {
            sum += list[i];
        }

        return sum / (double) list.length;
    }

    public static double getMedian(double[] values) {
        if (values.length == 1) {
            return select(values, values.length / 2 + 1);
        } else {
            double x = select(values, values.length / 2 + 1);
            double y = select(values, values.length / 2 + 1 + 1);

            return (x + y) / 2;
        }
    }

    public static double getMedian(ArrayList<Double> list) {
        if(list.size() == 0)
        {
            return Double.NaN;
        }
        
        if (list.size() == 1) {
            return list.get(0);
        }

        if (list.size() % 2 == 1) {
            return select(list, list.size() / 2 + 1);
        } else {
            Double x = select(list, list.size() / 2 + 1);
            Double y = select(list, list.size() / 2 + 1 + 1);

            return (x.doubleValue() + y.doubleValue()) / 2;
        }
    }
    
    

    public static double getPercentile(ArrayList<Double> list, double p) {
        if (list.size() == 1) {
            return list.get(0);
        } else {
            return select(list, (int) Math.round(p * (double) list.size()));
        }
    }

    public static double getMedian2(ArrayList<Double> list) {
        if (list.size() == 0) {
            return Double.NaN;
        }

        Collections.sort(list);
        if (list.size() % 2 == 1) {
            return list.get(list.size() / 2);
        } else {
            Double x = list.get(list.size() / 2);
            Double y = list.get(list.size() / 2 + 1);

            return (x.doubleValue() + y.doubleValue()) / 2;
        }
    }

    public static Double select(ArrayList<Double> list, int k) {
        if (list.size() <= 10) {
            Collections.sort(list);
            return list.get(Math.min(Math.max(0, k - 1), list.size() - 1));
        }

        int numSubsets = list.size() / 5;
        ArrayList[] subsets = new ArrayList[numSubsets];
        for (int i = 0; i < numSubsets; i++) {
            subsets[i] = new ArrayList<Double>();
            for (int j = 0; j < 5; j++) {
                subsets[i].add(list.get(i * 5 + j));
            }
        }

        ArrayList<Double> x = new ArrayList<Double>(numSubsets);
        for (int i = 0; i < numSubsets; i++) {
            x.add(select(subsets[i], 3));
        }


        Double M = select(x, list.size() / 10);

        ArrayList<Double> L1 = new ArrayList<Double>();
        ArrayList<Double> L2 = new ArrayList<Double>();
        ArrayList<Double> L3 = new ArrayList<Double>();
        for (int i = 0; i < list.size(); i++) {
            Double item = list.get(i);

            if (item.compareTo(M) < 0) {
                L1.add(item);
            } else if (item.compareTo(M) > 0) {
                L3.add(item);
            } else {
                L2.add(item);
            }
        }

        if (k <= L1.size()) {
            return select(L1, k);
        } else if (k > L1.size() + L2.size()) {
            return select(L3, k - L1.size() - L2.size());
        }

        return M;
    }
    
    public static double getMin(ArrayList<Double> list)
    {
    	double min = Double.MAX_VALUE;
    	for(double d : list)
    	{
    		min = Math.min(d, min);
    	}
    	return min;
    }
    
    public static double getMax(ArrayList<Double> list)
    {
    	double max = Double.MIN_VALUE;
    	for(double d : list)
    	{
    		max = Math.max(d, max);
    	}
    	return max;
    }

    public static double select(double[] list, int k) {
        if (list.length <= 10) {
            Arrays.sort(list);
            return list[Math.min(k - 1, list.length - 1)];
        }

        int numSubsets = list.length / 5;
        ArrayList[] subsets = new ArrayList[numSubsets];
        for (int i = 0; i < numSubsets; i++) {
            subsets[i] = new ArrayList<Double>();
            for (int j = 0; j < 5; j++) {
                subsets[i].add(list[i * 5 + j]);
            }
        }

        ArrayList<Double> x = new ArrayList<Double>(numSubsets);
        for (int i = 0; i < numSubsets; i++) {
            x.add(select(subsets[i], 3));
        }


        Double M = select(x, list.length / 10);

        ArrayList<Double> L1 = new ArrayList<Double>();
        ArrayList<Double> L2 = new ArrayList<Double>();
        ArrayList<Double> L3 = new ArrayList<Double>();
        for (int i = 0; i < list.length; i++) {
            Double item = list[i];

            if (item.compareTo(M) < 0) {
                L1.add(item);
            } else if (item.compareTo(M) > 0) {
                L3.add(item);
            } else {
                L2.add(item);
            }
        }

        if (k <= L1.size()) {
            return select(L1, k);
        } else if (k > L1.size() + L2.size()) {
            return select(L3, k - L1.size() - L2.size());
        }

        return M;
    }

    /**
     * Returns two-tailed p-value from z-value
     *
     * @param Z
     * @return
     */
    public static double NormalZ(double Z) {

        double Y, X, w, temp, Temp2;
        double Z_MAX, NormalZx = 0, WinP;

        Z_MAX = 6;

        if (Math.abs(Z) < 5.9999999) {
            if (Z == 0.0) {
                return 1;
            } else {
                Y = 0.5 * Math.abs(Z);
                if (Y >= (Z_MAX * 0.5)) {
                    X = 1.0;
                } else if (Y < 1.0) {

                    w = Y * Y;
                    X = ((((((((0.000124818987 * w - 0.001075204047) * w + 0.005198775019) * w - 0.019198292004) * w + 0.059054035642) * w - 0.151968751364) * w + 0.319152932694) * w - 0.5319230073) * w + 0.797884560593) * Y * 2.0;

                } else {

                    Y = Y - 2.0;
                    X = (((((((((((((-0.000045255659 * Y
                            + 0.00015252929) * Y - 0.000019538132) * Y
                            - 0.000676904986) * Y + 0.001390604284) * Y
                            - 0.00079462082) * Y - 0.002034254874) * Y
                            + 0.006549791214) * Y - 0.010557625006) * Y
                            + 0.011630447319) * Y - 0.009279453341) * Y
                            + 0.005353579108) * Y - 0.002141268741) * Y
                            + 0.000535310849) * Y + 0.999936657524;
                }

                if ((X + 1.0) < (1.0 - X)) {
                    NormalZx = (X + 1.0);
                } else {
                    NormalZx = (1.0 - X);
                }

            }
        } else {
            temp = ((Math.abs(Z) - 5.999999) * 10);
            Temp2 = Math.pow(1.6, temp);
            WinP = Math.pow(10, -9);
            WinP = WinP / Temp2;
            NormalZx = WinP;

        }
        return (NormalZx);
    }
}