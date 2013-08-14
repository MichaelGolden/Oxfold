package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class BasePairWeightingOxfoldNormalized implements BasePairWeighting{
	
	double[][] distances; 
	double weight;
	boolean [][] canPair;
	double normalization;
	
	public BasePairWeightingOxfoldNormalized(double [][] distances, double weight, boolean [][] canPair)
	{
		this.distances = distances;
		this.weight = weight;
		this.canPair = canPair;
		this.normalization = calculateNormalization();
	}
	
	public BasePairWeightingOxfoldNormalized(int [][] distances, double weight, boolean [][] canPair)
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
		this.canPair = canPair;
		this.normalization = calculateNormalization();
	}
	
	private double calculateNormalization()
	{
		double min = Double.MAX_VALUE;
		for(int i = 0 ; i < distances.length ; i++)
		{
			for(int j = i+1; j < distances.length ; j++)
			{
				if(canPair[i][j])
				{

					//min = Math.min(distances[i][j]/(weight*Math.abs(i-j)), min);
					min = Math.min(distances[i][j]/(weight), min);
				}
			}
		}
		
		return min;
	}

	@Override
	public double getDoubleWeight(int i, int j) {
		//System.out.println("Retrieving weight double");
		return Math.exp(-((distances[i][j]/(weight*Math.abs(i-j)))-normalization));
	}

	@Override
	public PointRes getPointResWeight(int i, int j) {
		//return PointRes.valueOf(Math.exp(-distances[i][j]/weight));
		//System.out.println("Retrieving weight PointRes");
		//double val = -((distances[i][j]/(weight*Math.abs(i-j)))-normalization);
		double val = -((distances[i][j]/(weight))-normalization);
		System.out.println("Calculating "+i+"\t"+j+"\t"+val+"\t"+Math.exp(val)+"\t"+normalization+"\t"+Math.exp(-distances[i][j]/0.5));
		return PointRes.valueOf(Math.exp(val));
	}
}
