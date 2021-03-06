package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;



public class Helix {
		private int leftIdx;
		private int rightIdx;
		private int helixLength;
		private PointRes score;
		public double maxdelta;
		
		public Helix(int leftIdx, int rightIdx, int helixLength) {
			this.leftIdx = leftIdx;
			this.rightIdx = rightIdx;
			this.helixLength = helixLength;
		}
		
		public Helix(int leftIdx, int rightIdx, int helixLength, PointRes[][] diffs) {
			this.leftIdx = leftIdx;
			this.rightIdx = rightIdx;
			this.helixLength = helixLength;
			PointRes score = PointRes.ZERO;
			for (int j = 0; j<helixLength; j++) {
				score = score.add(diffs[leftIdx+j][rightIdx-j]);
			}
			this.score = score;
		}
		
		public Helix(int leftIdx, int rightIdx, int helixLength, double[][] diffs) {
			this.leftIdx = leftIdx;
			this.rightIdx = rightIdx;
			this.helixLength = helixLength;
			double score = 0.0;
			for (int j = 0; j<helixLength; j++) {
				score = score+diffs[leftIdx+j][rightIdx-j];
			}
			this.score = PointRes.valueOf(score);
		}
		
		public Helix() {
			this.leftIdx = -1; this.rightIdx = -1;
			this.helixLength = 0; this.score = PointRes.ZERO;
		}
		
		public void setScore(PointRes score) {
			this.score = score;
		}
		
		public int getLeftIdx() {
			return leftIdx;
		}
		
		public int getRightIdx() {
			return rightIdx;
		}
		
		public int getHelixLength() {
			return helixLength;
		}
		
		public PointRes getScore() {
			return score;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + helixLength;
			result = prime * result + leftIdx;
			result = prime * result + rightIdx;
			result = prime * result + ((score == null) ? 0 : score.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Helix other = (Helix) obj;
			if (helixLength != other.helixLength)
				return false;
			if (leftIdx != other.leftIdx)
				return false;
			if (rightIdx != other.rightIdx)
				return false;
			return true;
		}

		@Override
		public String toString() {
			return "Helix [leftIdx=" + leftIdx + ", rightIdx=" + rightIdx
					+ ", helixLength=" + helixLength + "]";
		}
		
		

		
}
