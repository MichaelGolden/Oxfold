package uk.ac.ox.osscb.inoutside;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Arrays;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.PointRes;

public class PosteriorProbabilities {

		private PointRes[] unpairedProbs;
		private PointRes[][] pairedProbs;
		private PointRes maxP;
		private int maxLeftIdx;
		private int maxRightIdx;
		
		public PosteriorProbabilities(int [] structure)
		{
			pairedProbs = new PointRes[structure.length][structure.length];
			for(int i = 0 ; i < structure.length ; i++)
			{
				if(structure[i] != Constants.UnpairedBaseIdx)
				{
					pairedProbs[i][structure[i]] = PointRes.ONE;
				}
			}
		}
		
		public PosteriorProbabilities(PointRes[] unpairedProbs, PointRes[][] pairedProbs, PointRes maxP, int maxLeftIdx, int maxRightIdx) {
			this.unpairedProbs = unpairedProbs;
			this.pairedProbs = pairedProbs;
			this.maxP = maxP;
			this.maxLeftIdx = maxLeftIdx;
			this.maxRightIdx = maxRightIdx;
		}
				
		public PointRes getMaxP() {
			return maxP;
		}
		
		public PointRes[] getUnpairedProbs() {
			return unpairedProbs;
		}
		
		public PointRes[][] getPairedProbs() {
			return pairedProbs;
		}
		
		public int getMaxLeftIdx() {
			return maxLeftIdx;
		}
		
		public int getMaxRightIdx() {
			return maxRightIdx;
		}
		
		public double [][] getBasePairProbs()
		{
			double [][] basePairProbs = new double[pairedProbs.length][pairedProbs.length];
			for(int i = 0 ; i < pairedProbs.length ; i++)
			{
				for(int j = 0 ; j < pairedProbs.length ; j++)
				{
					if(pairedProbs[i][j].doubleValue() > 0)
					{
						basePairProbs[i][j]  = pairedProbs[i][j].doubleValue();
						basePairProbs[j][i] = basePairProbs[i][j];
					}
				}
			}
			return basePairProbs;
		}
		
		public void savePosteriorProbabilities(File outputFile) throws IOException
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
			
			/*
			PointRes [][] basePairProbs = new PointRes[pairedProbs.length][pairedProbs.length];
			PointRes [] totalProbs = new PointRes[pairedProbs.length];
			for(int i = 0 ; i < pairedProbs.length ; i++)
			{
				totalProbs[i] = new PointRes(1);
				for(int j = 0 ; j < pairedProbs.length ; j++)
				{
					totalProbs[i]=  totalProbs[i].subtract(pairedProbs[i][j]);
				}
				PointRes pairProb = PointRes.ONE.subtract(unpairedProbs[i]);
				totalProbs[i] = PointRes.ONE.subtract(totalProbs[i]);
				System.out.println("MMM "+totalProbs[i].doubleValue()+"\t"+pairProb.doubleValue());
			}*/
			
			for(int i = 0 ; i < pairedProbs.length ; i++)
			{
				for(int j = 0 ; j < pairedProbs[0].length ; j++)
				{
					if(pairedProbs[i][j] != null)
					{
						writer.write(pairedProbs[i][j].toString()+",");
					}
					else
					{
						writer.write("0.0,");
					}
				}
				writer.write("\n");
			}
			writer.close();
		}
}
