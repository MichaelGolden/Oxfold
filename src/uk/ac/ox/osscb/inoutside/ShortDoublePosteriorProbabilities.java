package uk.ac.ox.osscb.inoutside;


public class ShortDoublePosteriorProbabilities {
		
	private double[] unpairedProbs;
	private double[][] pairedProbs;
	
	public ShortDoublePosteriorProbabilities(double[] unpairedProbs, double[][] pairedProbs) {
		this.unpairedProbs = unpairedProbs;
		this.pairedProbs = pairedProbs;
	}
			
	public double[] getUnpairedProbs() {
		return unpairedProbs;
	}
	
	public double[][] getPairedProbs() {
		return pairedProbs;
	}
	
}
