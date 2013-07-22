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
		return alpha*(Math.exp(-Math.abs(i-j)/tau) - 1) + 1;
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		return PointRes.valueOf(alpha*(Math.exp(-Math.abs(i-j)/tau) - 1) + 1);
	}

}
