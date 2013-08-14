package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class BasePairWeightingCofold implements BasePairWeighting {
	
	double alpha;
	double tau;
	
	public BasePairWeightingCofold(double alpha, double tau)
	{
		this.alpha = alpha;
		this.tau = tau;
	}
	
	@Override
	public double getDoubleWeight(int i, int j) {
		//return alpha*(Math.exp(-Math.abs(i-j)/tau) - 1) + 1;

		return Math.exp(-Math.abs(i-j)/tau);
		//return Math.abs(i-j) > tau ? 0 : 1;
		
		//return Math.max(0.5, 1-Math.pow(Math.abs(i-j)/1000,4));
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		//return PointRes.valueOf(alpha*(Math.exp(-Math.abs(i-j)/tau) - 1) + 1);
		return PointRes.valueOf(Math.exp(-Math.abs(i-j)/tau));
		//return PointRes.valueOf(Math.max(0.5, 1-Math.pow(Math.abs(i-j)/1000,4)));
	}

}
