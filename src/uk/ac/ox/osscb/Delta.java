package uk.ac.ox.osscb;

import java.math.BigDecimal;

public class Delta {
	
	public static   Double IterationCutOffDouble;
	public   BigDecimal IterationCutOff = BigDecimal.valueOf(IterationCutOffDouble);
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
