package uk.ac.ox.osscb.inoutside;


public class PPOutputDouble {

	/**
	 * new Helix
	 */
	private Helix helix;
		
	/**
	 * Expected gain in predictive accuracy
	 */
	private double diff;
	
	public PPOutputDouble(Helix helix, double diff) {
		super();
		this.helix = helix;
		this.diff = diff;
	}

	public Helix getHelix() {
		return this.helix;
	}
	
	public int getLeftIdx() {
		return this.helix.getLeftIdx();
	}

	public int getRightIdx() {
		return this.helix.getRightIdx();
	}

	public double getDiff() {
		return diff;
	}
	
	public int gethelixLength() {
		return this.helix.getHelixLength();
	}

	
}
