package uk.ac.ox.osscb.inoutside;

import uk.ac.ox.osscb.PointRes;

public class PointResUpperMatrix {
	
	float [] fraction;
	int [] exponent;
	
	public int length;
	
	public PointResUpperMatrix(int length)
	{
		this.length = length;
		this.fraction = new float[(length*length+length)/2];
		this.exponent = new int[(length*length+length)/2];
	}
	
	public PointRes get(int i, int j)
	{
		int index = ((j*j+j)/2)+i;
		return new PointRes(fraction[index], exponent[index]);
	}
	
	public void set(int i, int j, PointRes value)
	{
		if(i > j)
		{
			return;
		}

		int index = ((j*j+j)/2)+i;
		this.fraction[index] = value.fraction;
		this.exponent[index] = value.exponent;
	}
}
