package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class BasePairWeightingOxfold implements BasePairWeighting{
	
	double[][] distances; 
	double weight;
	
	public BasePairWeightingOxfold(double [][] distances, double weight)
	{
		this.distances = distances;
		this.weight = weight;
	}
	
	public BasePairWeightingOxfold(int [][] distances, double weight)
	{
		this.distances = new double[distances.length][distances.length];
		for(int i = 0 ; i < distances.length ; i++)
		{
			for(int j = 0 ; j < distances.length ; j++)
			{
				this.distances[i][j] = distances[i][j];
			}
		}
		this.weight = weight;
	}

	@Override
	public double getDoubleWeight(int i, int j) {
		return Math.exp(-distances[i][j]/weight);
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		return PointRes.valueOf(Math.exp(-distances[i][j]/weight));
	}
}
