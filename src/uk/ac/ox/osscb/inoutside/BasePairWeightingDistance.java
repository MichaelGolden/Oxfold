package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class BasePairWeightingDistance implements BasePairWeighting {
	
	double tau;
	
	public BasePairWeightingDistance(double tau)
	{
		this.tau = tau;
	}
	
	@Override
	public double getDoubleWeight(int i, int j) {
		double alpha = 3.0;
		double beta = 0.5;
		
		//factor = Math.exp(-distance2/tau);
		return beta + (1-beta)*0.5*(-Math.tanh(alpha*((Math.abs(i-j)-tau)/tau)))+1;
		//return Math.exp(-Math.abs(i-j)/tau);
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		double alpha = 3.0;
		double beta = 0.5;
		return PointRes.valueOf(beta + (1-beta)*0.5*(-Math.tanh(alpha*((Math.abs(i-j)-tau)/tau)))+1);
		//return PointRes.valueOf(Math.exp(-Math.abs(i-j)/tau));
	}

}
