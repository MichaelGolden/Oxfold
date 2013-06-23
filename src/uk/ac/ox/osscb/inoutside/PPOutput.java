package uk.ac.ox.osscb.inoutside;

import java.math.BigDecimal;

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
	private BigDecimal diff;
	/**
	 * Compatibility score 
	 */
	private BigDecimal comp;
	
	/**
	 * Helix length
	 */
	private int helixLength;
	
	public PPOutput(int leftIdx, int rightIdx, int helixLength, BigDecimal diff, BigDecimal comp) {
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

	public BigDecimal getDiff() {
		return diff;
	}
	
	public BigDecimal getComp() {
		return comp;
	}
	
	public int gethelixLength() {
		return helixLength;
	}
}

