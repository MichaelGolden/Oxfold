package uk.ac.ox.osscb.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

import uk.ac.ox.osscb.DistancesCalculator;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.analysis.RNAFoldingTools.MultiThreadedPosteriorDecoding;
import uk.ac.ox.osscb.analysis.RNAFoldingTools.Pair;
import uk.ac.ox.osscb.analysis.RNAFoldingTools.MultiThreadedPosteriorDecoding.PosteriorDecodingThread;
import uk.ac.ox.osscb.inoutside.ParallelInsideOutsideCalculator.Sector;

public class HaasMEA {
	 /**
     * A single-threaded method for generating the posterior-decoding structure.
     *
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb an array of length N representing probabilities for
     * unpaired bases.
     * @return an array of paired positions. Where (i, array[i]) represents a
     * nucleotide pairing between nucleotides (i+1, array[i]), if array[i] = 0,
     * then (i+1) is unpaired.
     */
    public static int[] getPosteriorDecodingConsensusStructure(double[][] basePairProb, double[] singleBaseProb) {


        double[][] eMatrix = new double[basePairProb.length][basePairProb[0].length];
        for (int i = 0; i < eMatrix.length; i++) {
            for (int j = 0; j < eMatrix.length; j++) {
                eMatrix[i][j] = RNAFoldingTools.emptyValue;
            }
        }

        int[] pairedWith = new int[eMatrix.length];
        int[][] S = new int[basePairProb.length][basePairProb[0].length];
        recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, 0, eMatrix.length - 1);

        //printMatrix(eMatrix);
        //System.out.println();
        //printMatrix(S);
       // writeMatrix(S, new File("e.matrix"));
        //writeMatrix(S, new File("s.matrix"));
        traceBack(S, 0, eMatrix.length - 1, pairedWith);

        return pairedWith;
    }
    
    public static double [][] getDeltas(File bpFile) throws IOException
	{
		double [][] basePairProb = null;
		
		BufferedReader buffer = new BufferedReader(new FileReader(bpFile));
		String textline = null;
		int i = 0;
		while((textline = buffer.readLine()) != null)
		{
			String [] split = textline.split("[,(\\s)](\\s)*");
			if(i == 0)
			{
				basePairProb = new double[split.length][split.length];
			}
			
			for(int j = 0 ; j < basePairProb.length ; j++)
			{
				basePairProb[i][j] = Double.parseDouble(split[j]);
				if(basePairProb[i][j] != -1)
				{
					basePairProb[j][i] = basePairProb[i][j];
				}
			}
			
			i++;
		}
		
		for(int x = 0 ; x < basePairProb.length ; x++)
		{
			for(int y = 0 ; y < basePairProb.length ; y++)
			{
				if(basePairProb[x][y] != -1)
				{
					basePairProb[y][x] = basePairProb[x][y];
				}
			}
		}
		
		buffer.close();
		
		return basePairProb;
	}
    
    public static void print(double [][] mat)
    {
    	for(int i = 0 ; i < mat.length ;i++)
    	{
    		for(int j = 0 ; j < mat.length ;j++)
        	{
        		System.out.print(mat[i][j]+",");
        	}
    		System.out.println();
    	}
    }
    
    public static void print(int [][] mat)
    {
    	for(int i = 0 ; i < mat.length ;i++)
    	{
    		for(int j = 0 ; j < mat.length ;j++)
        	{
        		System.out.print(mat[i][j]+",");
        	}
    		System.out.println();
    	}
    }
    
    
    public static void main(String [] args)
    {
    	try {
			double [] [] bp = StructureData.getBasePairProb(new File("test.bp"));
			double [] single = StructureData.getSingleProbs(bp);
			double [] [] deltas = getDeltas(new File("test.deltas"));
			//print(bp);
			//System.out.println();
			//print(deltas);
			   
			double[][] eMatrix = new double[bp.length][bp[0].length];
	        for (int i = 0; i < eMatrix.length; i++) {
	            for (int j = 0; j < eMatrix.length; j++) {
	                eMatrix[i][j] = RNAFoldingTools.emptyValue;
	            }
	        }

	        int[] pairedWith = new int[eMatrix.length];
	        int[][] S = new int[eMatrix.length][eMatrix.length];    
	        for(int i = 0 ; i < S.length ; i++)
	        {
	        	Arrays.fill(S[i], -1);
	        }
	        
	    /*    for(int i = 0 ; i < eMatrix.length ; i++)
	        {
	        	for(int j = 0 ; j < eMatrix.length ; j++)
		        {
	        		haasRecurse(bp,single, deltas, eMatrix, S, i, j,0.0);
		        }	        	
	        }*/
	        double delta = 0.8;
	        double tau = Double.POSITIVE_INFINITY;
	      //  recursePosteriorDecoding(bp,single,eMatrix,S,0,eMatrix.length-1);
			haasRecurse(bp,single, deltas, eMatrix, S, 0, eMatrix.length-1,delta,tau);
			
	    
	       System.out.println();
	       print(eMatrix);
	       System.out.println("SMatrix");
	       print(S);

	        System.out.println("Delta="+delta);
	       //int [] pairedSites = new int[eMatrix.length];
	       try
	       {
		      traceBack(S, 0, eMatrix.length-1, pairedWith);
		      
		      System.out.println("A"+RNAFoldingTools.getDotBracketStringFromPairedSites(pairedWith));
		      int [] mea = RNAFoldingTools.getPosteriorDecodingConsensusStructure(bp);
		      System.out.println("B"+RNAFoldingTools.getDotBracketStringFromPairedSites(mea));
		      
		      for(int i = 0 ; i < mea.length ; i++)
		      {
		    	  if(mea[i] != 0)
		    	  {
		    		  int x = i;
		    		  int y = mea[i]-1;
		    		  System.out.println("D "+x+"\t"+y+"\t"+deltas[x][y]);
		    	  }
		      }

		      System.out.println("E.....<.<<<<<.<<<<....<<<<<..<<<...>>>>>>>>....>>>>>>>>>>...".replace('<', '(').replace('>', ')'));
	       }
	       catch(Exception ex)
	       {
	    	   System.err.println("Stackoverflow");
	       }
	    		 
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    
    public static final int CAN_PAIR = -2;

    /**
     * A recursive method that fills the dynamic programming matrix for
     * generating the posterior decoding structure.
     *
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb a vector of length N representing the probability
     * that a base at a position is unpaired.
     * @param eMatrix the dynamic programming matrix for the posterior-decoding.
     * @param i the start position of the window.
     * @param j the end position of the window.
     * @param pairedWith an array of paired positions. Where (i, array[i])
     * represents a pairing between nucleotides (i+1, array[i]), if array[i] =
     * 0, then (i+1) is unpaired.
     * @return the value of the eMatrix at position (i, j).
     */
    private static double haasRecurse2(double[][] basePairProb, double[] singleBaseProb, double [][] deltas, double[][] eMatrix, int[][] S, int i, int j, double delta) {    	
        if (i > j) {
        	System.out.println("i>j");
            return 0;
        }
        else
        if (eMatrix[i][j] != RNAFoldingTools.emptyValue) {
        	System.out.println("stored");
            return eMatrix[i][j];
        }
        else
        if (i == j) {
        	System.out.println("i==j");
            eMatrix[i][j] = 0.5*singleBaseProb[i];
            return eMatrix[i][j];
        }
        else // i < j
        {        	
        	double h6 = 0;
        	if(deltas[i][j] <= delta && S[i+1][j-1] == CAN_PAIR)
        	{
        		System.out.println("<delta and can pair");
	        	double h1 = Double.MIN_VALUE;
	        	int maxk = 0;
	        	for(int k = i ; k < j ; k++)
	        	{
	        		double bif = haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, i, k,delta)+haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, k+1, j,delta); 
	        		if(bif > h1)
	        		{
	        			h1 = bif;
	        			maxk = k;
	        		}        		
	        	}
	        	
	        	/*
	        	double pair =0;
	        	if(S[i+1][j-1] == CAN_PAIR)
	        	{
	        		pair = basePairProb[i][j]+ haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, i+1, j-1,delta);
	        	}*/
	        	double pair = basePairProb[i][j]+ haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, i+1, j-1,delta);

	        	if(pair > h1)
	        	{
	        		S[i][j] = CAN_PAIR;
	        		h1 = pair;
	        		System.out.println("A:"+i+","+j);
	        	}
	        	else
	        	{
	        		S[i][j] = maxk;
	        	}
	        	h6 = h1;
	        	
	        	eMatrix[i][j] = h6;
	        	return h6;
        	}
        	else
        	{        	
        		System.out.println(">delta");
	        	double h2 = Double.MIN_VALUE;
	        	int maxk = 0;
	        	for(int k = i ; k < j ; k++)
	        	{
	        		double bif = haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, i, k,delta)+haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, k+1, j,delta);
	        		if(bif > h2)
	        		{
	        			h2 = bif;
	        			maxk = k;
	        		}        		
	        	}
	        	
	        	double h3 = Double.MIN_VALUE;
	        	int maxkhelix = -1;
	        	if(deltas[i][j] > delta)
	        	{
		        	for(int k = 0 ; k <= (j-i)/2 - 1 ; k++)
		        	{
		        		double h4 = haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, i+k+1, j-k-1,delta);
		        		
		        		boolean indicator = false;
		        		
		        		for(int l = 0 ; l <= k ; l++) // does a strong base-pair exist in helix
		        		{
		        			if(deltas[i+l][k-l] > delta)
		        			{
		        				indicator = true;
		        				break;
		        			}
		        		}
		        		double h5 = 0;
		        		if(indicator)
		        		{
		        			for(int l = 0 ; l <= k ; l++)
		            		{
		        				h5 += basePairProb[i+l][j-l];		        			
		            		}
		        		}
		        		System.out.println("AAA"+k+"\t"+h4+"\t"+h5+"\t"+indicator);
		        		h4 = h4 + h5;
		        		if(h4 > h3)
		        		{
		        			h3 = h4;
		        			maxkhelix = k;
		        		}
		        	}
	        	}
	        	

        		System.out.println("DDD"+h3+"\t"+h2+"\t"+maxkhelix);
	        	if(h3 > h2)
	        	{
	        		//System.out.println(i+"\t"+j+"paired");
	        		//S[i][j] = CAN_PAIR;
	        		//for(int i )
	        		for(int l = 0 ; l <= maxkhelix ; l++)
	            	{
	        			S[i+l][j-l] = CAN_PAIR;
	            	}
	        	}
	        	else
	        	{
	        		S[i][j] = maxk;
	        		
	        	}
	        	
	        	h6 = Math.max(h2, h3);
	        	
	        	eMatrix[i][j] = h6;
	        	return h6;
        	}
        	//eMatrix[i][j] = h6;	
	        //return eMatrix[i][j];
        }
    }
    
    public static int [] haasMEA(PointRes[][] basePairProb, PointRes [] singleBaseProb, PointRes [][] deltas, double delta, double tau)
    {
    	double [][] bp = new double[basePairProb.length][basePairProb.length];
    	double [] singles = new double[basePairProb.length];
    	double [][] dts = new double[basePairProb.length][basePairProb.length];
    	
    	for(int i = 0 ; i < bp.length ; i++)
    	{
    		singles[i] = singleBaseProb[i].doubleValue();
    		for(int j = i ; j < bp.length ; j++)
        	{
        		bp[i][j] = basePairProb[i][j].doubleValue();
        		bp[j][i] = bp[i][j];
        		dts[i][j] =  deltas[i][j].doubleValue();
        		dts[j][i] = dts[i][j];
        	}
    	}
    	
    	return haasMEA(bp,singles,dts,delta,tau);
    }
    
    
    
    public static int [] haasMEA(double[][] basePairProb, double [] singleBaseProb, double [][] deltas, double delta, double tau)
    {
    	double[][] eMatrix = new double[basePairProb.length][basePairProb[0].length];
        for (int i = 0; i < eMatrix.length; i++) {
            for (int j = 0; j < eMatrix.length; j++) {
                eMatrix[i][j] = RNAFoldingTools.emptyValue;
            }
        }

        int[][] S = new int[eMatrix.length][eMatrix.length];    
        for(int i = 0 ; i < S.length ; i++)
        {
        	Arrays.fill(S[i], -1);
        }
        
        int [] storeSites = new int[S.length];
        for(int i = 0 ; i < storeSites.length ; i++)
        {
        	storeSites[i] = -1;
        }
        
       /* int [][] distance = new int[eMatrix.length][eMatrix.length];
        for (int i = 0; i < eMatrix.length; i++) {
            for (int j = 0; j < eMatrix.length; j++) {
            	distance[i][j] = Integer.MIN_VALUE;
            }
        }*/
        
		haasRecurse(basePairProb,singleBaseProb, deltas, eMatrix, S,  storeSites, 0, eMatrix.length-1,delta, tau);
		//print(S);
        int[] pairedWith = new int[eMatrix.length];
		traceBack(S, 0, eMatrix.length-1, pairedWith);
		
		return pairedWith;
    }
    
  /*  private static int [] decode(double[][] eMatrix, int i, int j)
    {
    	int [] pairedSites = new int[eMatrix.length];
    	
    }*/
    
    private static double haasRecurse(double[][] basePairProb, double[] singleBaseProb, double [][] deltas, double[][] eMatrix, int[][] S, int [] pairedSites, int i, int j, double delta, double tau) {    	
        if (i > j) {
            return 0;
        }
        else
        if (eMatrix[i][j] != RNAFoldingTools.emptyValue) {
            return eMatrix[i][j];
        }
        else
        if (i == j) {
            eMatrix[i][j] = 0.5*singleBaseProb[i];
            return eMatrix[i][j];
        }
        else // i < j
        {        	
        	double bifurcation = Double.MIN_VALUE;
        	int maxk = 0;
        	for(int k = i ; k < j ; k++)
        	{
        		double bif = haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, pairedSites, i, k,delta,tau)+haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, pairedSites, k+1, j,delta,tau); 
        		if(bif > bifurcation)
        		{
        			bifurcation = bif;
        			maxk = k;
        		}        		
        	}
        
        	double helix_formation = Double.MIN_VALUE;
        	int maxkhelix = -1;
        	for(int k = 0 ; k <= (j-i)/2 - 1 ; k++)
        	{
        		double h4 = haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, pairedSites, i+k+1, j-k-1,delta, tau);
        		
        		boolean indicator = false;
        		for(int l = 0 ; l <= k ; l++) // does a strong base-pair exist in helix
        		{
        			if(deltas[i+l][j-l] > delta)
        			{
        				indicator = true;
        				break;
        			}
        		}

        		if(indicator)
        		{
        			for(int l = 0 ; l <= k ; l++)
            		{        				
        				double factor = 1;
    					if(tau != Double.POSITIVE_INFINITY)
    					{
    						/*if(distance[i+l+1][j-1-l] == Integer.MIN_VALUE)
    						{
    							//distance[i+l+1][j-1-l] = DistancesCalculator.distCalc(pairedSites,i+l+1, j-1-l, distance);    
    						//	System.out.println("calculating d="+(i+l+1)+"\t"+(j-1-l));
    						}*/
    						/*else
    						{
    							if(distance[i+l+1][j-1-l] !=  DistancesCalculator.distCalc(pairedSites,i+l+1, j-1-l))
    							{
    								System.err.println((i+l+1)+"\t"+(j-1-l)+"\t"+"DD"+distance[i+l+1][j-1-l]+"\t"+ DistancesCalculator.distCalc(pairedSites,i+l+1, j-1-l));
    							}
    							else
    							{
    								System.out.println((i+l+1)+"\t"+(j-1-l)+"\t"+"DD"+distance[i+l+1][j-1-l]);
    							}
    						}*/
    						
        					//double d = DistancesCalculator.distCalc(pairedSites,i+l+1, j-1-l);
        					double alpha = 3.0;
        					double beta = 0.5;
        					factor = beta + (1-beta)*0.5*(-Math.tanh(alpha*(( DistancesCalculator.distCalc(pairedSites,i+l+1, j-1-l)-tau)/tau)))+1;
    					}
                		h4 += factor*basePairProb[i+l][j-l];
            		}
        			
            		if(h4 > helix_formation)
            		{
            			helix_formation = h4;
        				maxkhelix = k;
            		}
        		}
        	}        	
        	
        	if(helix_formation >= bifurcation && maxkhelix >= 0)
        	{
        		for(int k = 0 ; k <= maxkhelix ; k++)
            	{
        			S[i+k][j-k] = CAN_PAIR;
        			S[j-k][i+k] = CAN_PAIR;
        			pairedSites[i+k] = j-k;
            	}

            	eMatrix[i][j] = helix_formation;
        	}
        	else
        	{
        		S[i][j] = maxk;
        		eMatrix[i][j] = bifurcation;
        	}

	        return eMatrix[i][j];
        }
    }
    
    /**
     * A recursive method that fills the dynamic programming matrix for
     * generating the posterior decoding structure.
     *
     * @param basePairProb a NxN matrix of base-pairing probabilities.
     * @param singleBaseProb a vector of length N representing the probability
     * that a base at a position is unpaired.
     * @param eMatrix the dynamic programming matrix for the posterior-decoding.
     * @param i the start position of the window.
     * @param j the end position of the window.
     * @param pairedWith an array of paired positions. Where (i, array[i])
     * represents a pairing between nucleotides (i+1, array[i]), if array[i] =
     * 0, then (i+1) is unpaired.
     * @return the value of the eMatrix at position (i, j).
     */
    private static double recursePosteriorDecoding(double[][] basePairProb, double[] singleBaseProb, double[][] eMatrix, int[][] S, int i, int j) {
        if (i > j) {
            return 0;
        }

        if (eMatrix[i][j] != RNAFoldingTools.emptyValue) {
            return eMatrix[i][j];
        }

        if (i == j) {
            eMatrix[i][j] = singleBaseProb[i];
            return eMatrix[i][j];
        }


        double u1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i + 1, j) + singleBaseProb[i];
        double p1 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i + 1, j - 1) + 2 * basePairProb[i][j]; // * 2
        double u2 = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i, j - 1) + singleBaseProb[j];
        double p2 = 0;
        int max_k = i;
        for (int k = i; k < j; k++) {
            double p2k = recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, i, k) + recursePosteriorDecoding(basePairProb, singleBaseProb, eMatrix, S, k + 1, j);
            // remember K
            if (p2k > p2) {
                p2 = p2k;
                max_k = k;
            }
        }

        /*
         * double max = 0; if(p1 > u1 && p1 > p2 && p1 > u2) { max = p1; S[i][j]
         * = -3; } else if(u1 > p2 && u1 > u2) { max = u1; S[i][j] = -1; } else
         * if(p2 > u2) { max = p2; S[i][j] = max_k; } else { max = u2; S[i][j] =
         * -2; }
         */


        double max = 0;
        if (u1 > p1 && u1 > p2 && u1 > u2) {
            max = u1;
            S[i][j] = -1;
        } else if (u2 > p1 && u2 > p2) {
            max = u2;
            S[i][j] = -3;
        } else if (p1 > p2) {
            max = p1;
            S[i][j] = CAN_PAIR; // paired
        } else {
            max = p2;
            S[i][j] = max_k;
        }

        eMatrix[i][j] = max;

        return eMatrix[i][j];
    }
    
    public static void traceBack(int[][] S, int i, int j, int[] pairedWith) {
    	//System.out.println(i+"\t"+j+"\t"+S[i][j]);
        if (i >= j) {
            // do nothing
        } else if (S[i][j] == CAN_PAIR) {
            pairedWith[i] = j + 1;
            pairedWith[j] = i + 1;
            traceBack(S, i + 1, j - 1, pairedWith);
        } else {
            traceBack(S, i, S[i][j], pairedWith);
            traceBack(S, S[i][j] + 1, j, pairedWith);
        }
    }
    
    
    public static int [] parallelHaasMEA(PointRes[][] basePairProb, PointRes[] singleBaseProb, PointRes [][] deltas, double delta, double tau)
    {
    	double [][] bp = new double[basePairProb.length][basePairProb.length];
    	double [] singles = new double[basePairProb.length];
    	double [][] dts = new double[basePairProb.length][basePairProb.length];
    	
    	for(int i = 0 ; i < bp.length ; i++)
    	{
    		singles[i] = singleBaseProb[i].doubleValue();
    		for(int j = i ; j < bp.length ; j++)
        	{
        		bp[i][j] = basePairProb[i][j].doubleValue();
        		bp[j][i] = bp[i][j];
        		dts[i][j] =  deltas[i][j].doubleValue();
        		dts[j][i] = dts[i][j];
        	}
    	}
    	
    	
    	ParallelHaasMEA mea = new ParallelHaasMEA(bp, singles, dts, delta, tau);
    	mea.execute();
    	return mea.pairedWith;
    }
    
   
    
    public static class ParallelHaasMEA {

        double[][] basePairProb;
        double[] singleBaseProb;
        int[] pairedWith;
        double delta;
        double tau;
        final Integer lock = new Integer(0);
        int maxThreads = Runtime.getRuntime().availableProcessors();
        int divisions = maxThreads*3;
        int threadsUsed = 0;
        int currentBlockY = 0;
        int currentX = 0;
        int currentY = 0;
        int length;
        int startingDivisionSize;
        int divisionSize;
        int blockIncrement;
        double[][] eMatrix;
        double[][] deltas;
        int[][] S;
        int [] storeSites;
        RunInfo runInfo = new RunInfo();

        public ParallelHaasMEA(double[][] basePairCount, double[] singleBaseCount, double [][] deltas, double delta, double tau) {
            this.length = singleBaseCount.length;
            this.pairedWith = new int[this.length];
            this.basePairProb = basePairCount;
            this.singleBaseProb = singleBaseCount;
            this.divisionSize = Math.max(this.length / divisions, 10);
            this.blockIncrement = this.divisionSize;

            this.startingDivisionSize = this.divisionSize;
            //System.out.println(divisionSize);
            this.eMatrix = new double[length][length];
            for (int i = 0; i < eMatrix.length; i++) {
                for (int j = 0; j < eMatrix.length; j++) {
                    eMatrix[i][j] = RNAFoldingTools.emptyValue;
                }
            }
            this.deltas = deltas;
            this.delta = delta;
            this.tau = tau;
            this.S = new int[length][length];
            for(int i = 0 ; i < S.length ; i++)
            {
            	Arrays.fill(S[i], -1);
            }
            
            storeSites = new int[S.length];
            for(int i = 0 ; i < storeSites.length ; i++)
            {
            	storeSites[i] = -1;
            }
        }
        
        /**
         * Generates a sequence of coordinates of windows on which
         * posterior-decoding can be performed independently in parallel.
         *
         * @return the positions of the next window to perform
         * posterior-decoding.
         */
        public Pair getNextSection() {
            Pair p = getNextSectionRecursive();
            if (p.x == -1) {
                return p;
            } else {
                return new Pair(Math.min(p.x, this.length - 1), Math.min(p.y, this.length - 1));
            }
        }

        /**
         * @see ParallelHaasMEA#getNextSection()
         */
        private Pair getNextSectionRecursive() {
            int newX = this.currentX;
            int newY = this.currentY;

            synchronized (lock) {
                if (newX == -1 && newY == -1) {
                    return new Pair(-1, -1);
                } else if (newX == 0 && newY == 0) {
                    newX = 0;
                    newY = divisionSize;
                } else {
                    newX += divisionSize;
                    newY += divisionSize;


                    if (currentBlockY >= this.length) {
                        newX = -1;
                        newY = -1;
                    } else if (newX >= this.length) {
                        currentBlockY += this.blockIncrement;

                        divisionSize = (int) Math.max((double) (length - currentBlockY) / (double) maxThreads, 10);
        
                        newX = 0;
                        newY = currentBlockY + divisionSize;
                    }
                }
            }

            this.currentX = newX;
            this.currentY = newY;

            if (newY - divisionSize > length) {
                return getNextSection();
            }

            return new Pair(newX, newY);
        }
        
        public void execute()
        {
        	try {
				process();
				
				// do final recursion
				HaasMEA.haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, pairedWith, 0, eMatrix.length-1, delta, tau);
				
				pairedWith = new int[eMatrix.length];
				traceBack(S, 0, eMatrix.length-1, pairedWith);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
        
        public class ParallelInfo
    	{
    		private int upto = 0;
    		List<Pair> sectors;

    		public ParallelInfo(List<Pair> sectors)
    		{
    			this.sectors = sectors;
    		}

    		public Pair getMostCompleteSector()
    		{
    			for(int i = upto+1 ; i < sectors.size() ; i++)
    			{
    				if(sectors.get(i).finished)
    				{
    					upto = i;
    				}
    				else
    				{
    					break;
    				}
    			}

    			return sectors.get(upto);
    		}


    		public boolean canRun(Pair s)
    		{
    			Pair upto = getMostCompleteSector();
    		//	System.out.println("GE"+(upto.n+maxThreads)+"\t"+s.n);
    			if(upto.n +divisions < s.n)
    			{
    				//System.out.println("here");
    				return false;
    			}
    			return true;
    		}
    	}
        
        public List<Object> process()
    			throws InterruptedException, ExecutionException {

    		ExecutorService service = Executors.newFixedThreadPool(maxThreads);

    		List<Future<Object>> futures = new ArrayList<Future<Object>>();
    		final ArrayList<Pair> jobs = new ArrayList<Pair>();
    		for(int n = 0 ; ; n++)
    		{
                final Pair p = getNextSection();
                if(p.x == -1 && p.y == -1)
                {
                	break;
                }

            	p.n = n;
                jobs.add(p);
    		}
    		System.out.println(jobs);
    		
    		final ParallelInfo info = new ParallelInfo(jobs);
    		//final long startTime = System.nanoTime();
    		for(final Pair p : jobs)
    		{
    			Callable<Object> callable = new Callable<Object>() {
    				public Object call() throws Exception {
    					
    					while(!info.canRun(p))
    					{
    						
    					}
    					
    				//	System.out.println("Running "+p);
    					haasRecurse(basePairProb, singleBaseProb, deltas, eMatrix, S, storeSites, p.x, p.y, delta, tau);
    					//long endTime = System.nanoTime();
    					//double elapsed = (endTime-startTime)/1e9;
    					//System.out.println("Elapsed:"+elapsed);
    					p.finished = true;
    					return null;
    				}
    			};
    			futures.add(service.submit(callable));
    		}
    		   		

    		service.shutdown();


    		List<Object> outputs = new ArrayList<Object>();
    		for (Future<Object> future : futures) {
    			outputs.add(future.get());
    		}
    		return outputs;
    	}

        public int[] getPairedSites() {
            return pairedWith;
        }
    }

}
