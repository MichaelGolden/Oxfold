package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;



/**
 * Sketch. To be changed later
 * @author Vladimir, Pierre
 *
 */
public class PPOutput {
	
	/**
	 * Left index of the new base pair
	 */
	private int leftIdx;
		
	/**
	 * Right index of the new base pair
	 */
	private int rightIdx;
		
	/**
	 * Expected gain in predictive accuracy
	 */
	private PointRes diff;
	/**
	 * Compatibility score 
	 */
	private PointRes comp;
	
	/**
	 * Helix length
	 */
	private int helixLength;
	
	public PPOutput(int leftIdx, int rightIdx, int helixLength, PointRes diff, PointRes comp) {
		super();
		this.leftIdx = leftIdx;
		this.rightIdx = rightIdx;
		this.helixLength = helixLength;
		this.diff = diff;
		this.comp = comp;
	}

	public int getLeftIdx() {
		return leftIdx;
	}

	public int getRightIdx() {
		return rightIdx;
	}

	public PointRes getDiff() {
		return diff;
	}
	
	public PointRes getComp() {
		return comp;
	}
	
	public int gethelixLength() {
		return helixLength;
	}

	@Override
	public String toString() {
		return "PPOutput [leftIdx=" + leftIdx + ", rightIdx=" + rightIdx
				+ ", diff=" + diff + ", comp=" + comp + ", helixLength="
				+ helixLength + "]";
	}
	
	
}

