package uk.ac.ox.osscb;



public class Delta {
	
	public static   Double IterationCutOffDouble;
	public   PointRes IterationCutOff = PointRes.valueOf(IterationCutOffDouble);
	public Delta(){
		IterationCutOffDouble = 0.5; 
	}
	public Double getIterationCutOffDouble(){
		return IterationCutOffDouble.doubleValue();
	}
	
	public void setIterationCutOffDouble(double d){
		IterationCutOffDouble = d;
	}
}
