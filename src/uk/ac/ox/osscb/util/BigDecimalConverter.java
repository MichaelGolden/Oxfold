package uk.ac.ox.osscb.util;

import java.math.BigDecimal;
import java.math.BigInteger;

public class BigDecimalConverter {
	
	public static BigDecimal valueOf(double dbl, int digits){
		
		if(0 == dbl)
			return BigDecimal.ZERO;
		if(1.0 == dbl)
			return BigDecimal.ONE;
		if(10.0 == dbl)
			return BigDecimal.TEN;
		
		if(dbl > 0)
			return valueOfPositive(dbl, digits);
		
		return valueOfPositive(-dbl, digits).negate();
	}
	
	private static BigDecimal valueOfPositive(double dbl, int digits){
		if(dbl <= 0){
			throw new IllegalArgumentException(String.format(
					"dbl must be positive. Got: %g", dbl));
		}
		int scale = (int)Math.floor(Math.log10(dbl));
		int scaleShifted = scale - digits + 1;
		double dblScaled = dbl / Math.pow(10, scaleShifted);
		long dblRounded = Math.round(dblScaled);// cut off unnecessary digits
		
		BigInteger valueOf = BigInteger.valueOf(dblRounded);
		BigDecimal bigDecimal = new BigDecimal(valueOf, -scaleShifted);
		return bigDecimal;
	}
}
