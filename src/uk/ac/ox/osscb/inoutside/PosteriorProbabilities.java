package uk.ac.ox.osscb.inoutside;

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
}
