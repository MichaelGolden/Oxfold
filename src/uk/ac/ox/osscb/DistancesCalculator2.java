package uk.ac.ox.osscb;

import java.math.BigDecimal;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ShortDoublePosteriorProbabilities;

public class DistancesCalculator2 {
	/**
	 * calculates distances between bases in alignment
	 * @author lepuslapis, lmath
	 * @param structure of alignment: if s[j]=k then j is paired with k; unpaired bases have s[i]=-1;
	 * @link {@link Constants#UnpairedBaseIdx} 
	 * @return distances between bases in alignment
	 */
	public double[][] distCalc(int[] structure, PosteriorProbabilities ppProbs) {
		if (null==structure) {
			throw new IllegalArgumentException("Structure cannot be null.");
		}

		int length = structure.length;
		double[][] distances = new double[length][length];
		double[][] pairedProbsSums = getProbPairsOutsideMatrix(length, ppProbs);
		
		for (int b=1; b<length; b++) {
			for (int j=0; j<length-b; j++) {
				if ((structure[j]>j)&&(structure[j]<=j+b)) {
					distances[j][j+b]=distances[structure[j]][j+b];
				} else {
					distances[j][j+b]= pairedProbsSums[j][j+b]+distances[j+1][j+b];
				}
				
			}
		}
		
		
		for (int b=1; b<length; b++) 
			for(int j=0; j<length-b; j++)
				distances[j][j+b] = Math.abs(distances[j][j+b]/(double) b);
	
		
		return distances;
	}
	
	
	public double[][] distCalc(int[] structure, ShortDoublePosteriorProbabilities ppProbs) {
		if (null==structure) {
			throw new IllegalArgumentException("Structure cannot be null.");
		}

		int length = structure.length;
		double[][] distances = new double[length][length];
		double[][] pairedProbsSums = getProbPairsOutsideMatrix(length, ppProbs);
		
		for (int b=1; b<length; b++) {
			for (int j=0; j<length-b; j++) {
				if ((structure[j]>j)&&(structure[j]<=j+b)) {
					distances[j][j+b]=distances[structure[j]][j+b];
				} else {
					distances[j][j+b]= pairedProbsSums[j][j+b]+distances[j+1][j+b];
				}
				
			}
		}
		
		
		for (int b=1; b<length; b++) 
			for(int j=0; j<length-b; j++)
				distances[j][j+b] = Math.abs(distances[j][j+b]/(double) b);
	
		
		return distances;
	}
	
/*
	public double[][][] getProbPairsOutsideMatrix(int length, PosteriorProbabilities ppProbs){
	
		double[][][] pairOutsideProbs = new double[length][length][length];
		
		for(int b=1; b<length; b++)
			for(int i=1; i<length-b-1; i++)
				for(int k=0; k<length; k++){
					if(b<=2)
						pairOutsideProbs[k][i][i+b] = 1-ppProbs.getUnpairedProbs()[k].doubleValue();
					else{		
						double iPair, jPair;
						if(i+1<k)
							iPair = ppProbs.getPairedProbs()[i+1][k].doubleValue(); 
						else
							iPair = ppProbs.getPairedProbs()[k][i+1].doubleValue();
						
						if(k<i+b-1)
							jPair = ppProbs.getPairedProbs()[k][i+b-1].doubleValue();
						else
							jPair = ppProbs.getPairedProbs()[i+b-1][k].doubleValue();
						
						pairOutsideProbs[k][i][i+b] = pairOutsideProbs[k][i+1][i+b-1] - iPair - jPair;
					}
	
				}


		return pairOutsideProbs;
	}*/
	
	public double[][] getProbPairsOutsideMatrix(int length, PosteriorProbabilities ppProbs){
		
		double[][] pairOutsideProbs = new double[length][length];
		
		for(int i=1; i<length-1; i++) {
			pairOutsideProbs[i][i+1] = 1-ppProbs.getUnpairedProbs()[i+1].doubleValue(); 
		}
		for(int b=2; b<length; b++) {
			for(int i=1; i<length-b; i++) {
				pairOutsideProbs[i][i+b] = pairOutsideProbs[i][i+b-1] - ppProbs.getPairedProbs()[i+1][i+b-1].doubleValue(); 
			}	
		}
			
		return pairOutsideProbs;
	}
	
	public double[][] getProbPairsOutsideMatrix(int length, ShortDoublePosteriorProbabilities ppProbs){
		
		double[][] pairOutsideProbs = new double[length][length];
		
		for(int i=1; i<length-1; i++) {
			pairOutsideProbs[i][i+1] = 1-ppProbs.getUnpairedProbs()[i+1]; 
		}
		for(int b=2; b<length; b++) {
			for(int i=1; i<length-b; i++) {
				pairOutsideProbs[i][i+b] = pairOutsideProbs[i][i+b-1] - ppProbs.getPairedProbs()[i+1][i+b-1]; 
			}	
		}
			
		return pairOutsideProbs;
	}
	
}
