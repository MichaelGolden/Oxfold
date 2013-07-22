package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public interface BasePairWeighting {
	public double getDoubleWeight(int i , int j);
	public PointRes getPointResWeight(int i , int j);
}
