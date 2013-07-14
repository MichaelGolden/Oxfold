package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;



public class PPOutputHelix {
	
		private Helix helix;
		private PointRes diff;
		private PointRes[][] diffs;
		private PointRes comp;
		
		public PPOutputHelix(Helix helix, PointRes[][] diffs, PointRes diff, PointRes comp) {
			this.helix = helix;
			this.diffs = diffs;
			this.diff = diff;
			this.comp = comp;
		}
		
		public Helix getHelix() {
			return this.helix;
		}
		
		public PointRes getDiff() {
			return this.diff;
		}
		
		public PointRes[][] getDiffs() {
			return this.diffs;
		}
		
		public PointRes getComp() {
			return this.comp;
		}
		
}
