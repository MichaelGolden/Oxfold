package uk.ac.ox.osscb.inoutside;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;

public class PosteriorProbabilities {

		private BigDecimal[] unpairedProbs;
		private BigDecimal[][] pairedProbs;
		private BigDecimal maxP;
		private int maxLeftIdx;
		private int maxRightIdx;
		
		public PosteriorProbabilities(BigDecimal[] unpairedProbs, BigDecimal[][] pairedProbs, BigDecimal maxP, int maxLeftIdx, int maxRightIdx) {
			this.unpairedProbs = unpairedProbs;
			this.pairedProbs = pairedProbs;
			this.maxP = maxP;
			this.maxLeftIdx = maxLeftIdx;
			this.maxRightIdx = maxRightIdx;
		}
				
		public BigDecimal getMaxP() {
			return maxP;
		}
		
		public BigDecimal[] getUnpairedProbs() {
			return unpairedProbs;
		}
		
		public BigDecimal[][] getPairedProbs() {
			return pairedProbs;
		}
		
		public int getMaxLeftIdx() {
			return maxLeftIdx;
		}
		
		public int getMaxRightIdx() {
			return maxRightIdx;
		}
		
		public void savePosteriorProbabilities(File outputFile) throws IOException
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
			for(int i = 0 ; i < pairedProbs.length ; i++)
			{
				for(int j = 0 ; j < pairedProbs[0].length ; j++)
				{
					writer.write(pairedProbs[i][j].doubleValue()+", ");
				}
				writer.write("\n");
			}
			writer.close();
		}
}
