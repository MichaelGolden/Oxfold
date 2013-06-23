package uk.ac.ox.osscb.inoutside;

import java.math.BigDecimal;

public class PPOutputHelix {
	
		private Helix helix;
		private BigDecimal diff;
		private BigDecimal[][] diffs;
		private BigDecimal comp;
		
		public PPOutputHelix(Helix helix, BigDecimal[][] diffs, BigDecimal diff, BigDecimal comp) {
			this.helix = helix;
			this.diffs = diffs;
			this.diff = diff;
			this.comp = comp;
		}
		
		public Helix getHelix() {
			return this.helix;
		}
		
		public BigDecimal getDiff() {
			return this.diff;
		}
		
		public BigDecimal[][] getDiffs() {
			return this.diffs;
		}
		
		public BigDecimal getComp() {
			return this.comp;
		}
		
}
