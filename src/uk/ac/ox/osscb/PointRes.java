package uk.ac.ox.osscb;



import java.io.Serializable;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * Extended exponent datatype definition
 * 
 * @author Z.Sukosd
 */

public class PointRes  extends Number implements Serializable {
	
	public static final PointRes ZERO = new PointRes(0, 0);
	public static final PointRes ONE = new PointRes(1);

	private static final long serialVersionUID = 6121110456253627605L;
	private static final double LOG_TWO = Math.log(2);
	// container for new number representation
	public float fraction; // = 0 OR >=1 0R <basis (in absolute value)
	public int exponent; // can be any integer, positive and negative

	public PointRes(float fraction, int exponent) {
		// Input fraction should always have exponent zero!!
		this.fraction = fraction;
		this.exponent = exponent;
	}
	
	public PointRes(PointRes point) {
		// Input fraction should always have exponent zero!!
		this.fraction = point.fraction;
		this.exponent = point.exponent;
	}

	public void setToFloat(float number) {
		this.exponent = 0;
		this.fraction = number;
		this.convert();
	}
	
	public void setToDouble(double number){
		long bits = Double.doubleToRawLongBits(number);
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		this.fraction = (float) fract;
		if (this.fraction == 0f) {
			this.fraction = 0;
			this.exponent = 0;
		} else {
			this.exponent = (int) ((bits >> 52) & 0x7FF) - 1023;
		}
	}

	public PointRes(double val) {
		long bits = Double.doubleToRawLongBits(val);
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		this.fraction = (float) fract;
		if (this.fraction == 0f) {
			this.fraction = 0;
			this.exponent = 0;
		} else {
			this.exponent = (int) ((bits >> 52) & 0x7FF) - 1023;
		}
	}
	
	public static PointRes max(PointRes a, PointRes b)
	{
		PointRes max = a;
		if(max.compareTo(b) < 0)
		{
			max = b;
		}
		return max;
	}

	@Override
	public PointRes clone() {
		PointRes newpoint = new PointRes(this.fraction, this.exponent);
		return newpoint.convert();
	}

	public void copyFrom(PointRes point) {
		this.exponent = point.exponent;
		this.fraction = point.fraction;
	}

	public int getExponent() {
		return exponent;
	}

	public void setExponent(int exponent) {
		this.exponent = exponent;
	}

	public float getFraction() {
		return fraction;
	}

	public void setFraction(float val) {
		this.fraction = val;
	}

	public PointRes convert() {
		if (fraction == 0f) {
			return this;
		}
		int bits = Float.floatToRawIntBits(fraction);
		int mantissaexp = ((bits >> 23) & 0xFF) - 127;
		exponent += mantissaexp;
		// set exp = 1 and then convert to double:
		fraction = Float.intBitsToFloat((bits & 0xBFFFFFFF) | 0x3f800000);
		return this;
	}

	public PointRes takeLogE(){
		double fraclog = Math.log(this.fraction);
		double exponentlog = LOG_TWO*this.exponent;
		this.fraction = (float) (fraclog + exponentlog);
		this.exponent = 0;
		this.convert();
		return this; 
	}
	
	public PointRes takeLog2(){
		double fraclog = Math.log(this.fraction)/LOG_TWO;
		this.fraction = (float) (fraclog + this.exponent);
		this.exponent = 0;
		this.convert();
		return this; 
	}
	
	public PointRes add(double d) {
		PointRes ret = this.clone();
		ret.fraction = ret.fraction + (float) d;
		ret.convert();
		return ret;
	}

	public PointRes add(PointRes point) {
		PointRes ret = this.clone();
		int diff = ret.exponent - point.exponent;
		if ((diff < -126 && point.fraction != 0) || ret.fraction == 0) {
			ret.exponent = point.exponent;
			ret.fraction = point.fraction;
		} else if ((diff > 127 && ret.fraction != 0) || point.fraction == 0) {
			// do nothing, adding zero to the number
		} else if (point.exponent != ret.exponent) {
			float newmantissa = point.fraction
					+ setExponent(ret.fraction, diff);
			ret.exponent = point.exponent;
			ret.fraction = newmantissa;
			ret.convert();
		} else {// equal exponents
			ret.fraction = ret.fraction + point.fraction;
			ret.convert();
		}
		return ret;
	}

	public boolean isSignificantlyLessThanZero(){
		if(this.fraction < 0f && this.exponent > -8){
			return true;
		}
		else return false; 
	}
	
	public int compareTo(PointRes point)
	{
		if(this.equals(point))
		{
			return 0;
		}
		else
		if(this.isLessThan(point))
		{
			return -1;
		}
		else
		{
			return 1;
		}
	}
	
	public boolean isLessThan(PointRes point) {
		if (this.fraction > 0f && point.fraction > 0f || this.fraction < 0f
				&& point.fraction < 0f) {
			if (this.exponent > point.exponent) {
				return false;
			} else if (this.exponent < point.exponent) {
				return true;
			} else { // case of equal exponents
				if (this.fraction < point.fraction) {
					return true;
				} else {
					return false;
				}
			}
		} else if (this.fraction >= 0f && point.fraction <= 0f) {
			return false;
		} else {
			return true;
		}
	}

	public boolean equals(PointRes point) {
		if (this.exponent == point.exponent && this.fraction == point.fraction) {
			//System.out.println(this.exponent +"\t"+point.exponent+"\t"+this.fraction+"\t"+point.fraction+" >>>1");
			return true;
		} else if (this.fraction == 0 && point.fraction == 0) {
			//System.out.println(point+" >>>2");
			return true;
		} else {
			//System.out.println(point+" >>>3");
			return false;
		}
	}

	public PointRes multiply(PointRes point) {
		PointRes ret = this.clone();
		ret.fraction = ret.fraction * point.fraction;
		ret.exponent = ret.exponent + point.exponent;
		ret.convert();
		return ret;
	}

	public PointRes multiply(double prob) {
		PointRes ret = this.clone();
		
		long bits = Double.doubleToRawLongBits(prob);
		// extract fraction part of probability
		double fract = Double
				.longBitsToDouble((bits & 0xBfffffffffffffffL) | 0x3ff0000000000000L);
		// cast to float
		float tmpfraction = (float) fract;
		if (tmpfraction == 0f) {
			// multiplying by zero
			ret.fraction = 0;
			ret.exponent = 0;
		} else {
			// multiplying by a non-zero value
			int tmpexponent = (int) ((bits >> 52) & 0x7FF) - 1023;
			ret.fraction = ret.fraction * tmpfraction;
			ret.exponent = ret.exponent + tmpexponent;
			ret.convert();
		}
		return ret;
	}

	public PointRes multiply(PointRes point, double prob) {
		PointRes ret = this.clone();
		ret.multiply(point);
		ret.multiply(prob);
		return ret;
	}

	public PointRes divide(PointRes point) {
		if (point.fraction == 0f) {
			System.err.println("Division by zero!");
			return null;
		}
		PointRes ret = this.clone();
		
		ret.fraction = ret.fraction / point.fraction;
		ret.exponent = ret.exponent - point.exponent;
		ret.convert();
		//System.out.println("divide"+ret);
		return ret;
	}
	
	public PointRes divide(PointRes point, int roundingMode) {
		return this.divide(point).round(new MathContext(512, RoundingMode.valueOf(roundingMode)));
	}
	
	public PointRes divide(PointRes point, RoundingMode roundingMode) {
		return this.divide(point).round(new MathContext(512, roundingMode));
	}
	
	public PointRes round(MathContext mc)
	{
		//System.out.println("round"+this.doubleValue());
		//BigDecimal bd = new BigDecimal(2);
		//System.out.println("A"+bd);
		//System.out.println("B"+this.exponent);
		//System.out.println("C"+this.fraction);
		//bd = bd.pow(this.exponent,mc).multiply(BigDecimal.valueOf(this.fraction));
		//bd = bd.round(mc);
		//return new PointRes(bd.doubleValue());
		return this.clone();
	}

	public float toFloat() {
		if (fraction == 0f) {
			return 0;
		} else if (-126 <= exponent && exponent <= 127) {
			// normal value
			return setExponent(fraction, exponent);
		} else if (exponent < -126) {
			return 0;
		} else {
			return fraction > 0f ? Float.POSITIVE_INFINITY
					: Float.NEGATIVE_INFINITY;
		}

	}

	public double toDouble() {
		if (fraction == 0f) {
			return 0;
		} else if (-126 <= exponent && exponent <= 127) {
			// normal value
			return (double) setExponent(fraction, exponent);
		} else if (exponent < -126) {
			return 0;
		} else {
			return fraction > 0f ? Double.POSITIVE_INFINITY
					: Double.NEGATIVE_INFINITY;
		}
	}

	public String toString() {
		return this.doubleValue()+"";
		//
	}
	
	public String getStringRepresentation()
	{
		return "" + this.fraction + " x 2^" +this.exponent;
	}

	// exponent must be between -126 and 127
	private static float setExponent(float d, int expin) {
		int bits = Float.floatToRawIntBits(d);
		bits = bits & 0x807fffff; // exp bits are 0
		bits = bits | ((expin + 127) << 23);
		return Float.intBitsToFloat(bits);
	}

	public void subtract(PointRes pr, PointRes tmp) {
		tmp.copyFrom(pr);
		tmp.multiply(-1);
		this.add(tmp);
	}
	
	public PointRes subtract(PointRes subtrahend)
	{
		PointRes a = this.clone();
		PointRes b = subtrahend.clone().multiply(-1);
				
		return a.add(b);
	}
	
	public static PointRes valueOf(float value)
	{
		return new PointRes(value);
	}
	
	public static PointRes valueOf(double value)
	{
		return new PointRes(value);
	}
	
	public int signum()
	{
		if(this.equals(PointRes.ZERO))
		{
			return 0;
		}
		else
		if(this.isLessThan(PointRes.ZERO))
		{
			return -1;
		}
		
		else
		{
			return 1;
		}
	}
	
	public static void main(String [] args)
	{
		PointRes a = new PointRes(10.5);
		PointRes b = new PointRes(4.3);
		System.out.println(a.multiply(b).toDouble());
		System.out.println(a.toDouble());
		System.out.println(b.toDouble());
		
		PointRes x = new PointRes(10.5);
		PointRes y = new PointRes(4.3);
		System.out.println(x.add(y));
		System.out.println(x);
		System.out.println(y);
		
	}

	@Override
	public double doubleValue() {		
		return this.toDouble();
	}

	@Override
	public float floatValue() {
		// TODO Auto-generated method stub
		return this.toFloat();
	}

	@Override
	public int intValue() {
		// TODO Auto-generated method stub
		return (int)this.toDouble();
	}

	@Override
	public long longValue() {
		// TODO Auto-generated method stub
		return (long)this.toDouble();
	}
}