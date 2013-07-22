package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class BasePairWeightingNull implements BasePairWeighting {

	@Override
	public double getDoubleWeight(int i, int j) {		
		return 1;
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		return PointRes.ONE;
	}

}
