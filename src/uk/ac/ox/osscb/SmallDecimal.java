package uk.ac.ox.osscb;

public class SmallDecimal {
	private double base;
	private long scale;
	
	public SmallDecimal(double val) {
		this.scale = (long) Math.floor(Math.log(Math.abs(val))/Math.log(2));
		this.base = val * Math.pow(2,-scale);
	}
	
	public SmallDecimal(double b, long s) {
		this.base = b; this.scale = s;
	}
	
	public long getScale() {
		return this.scale;
	}
	
	public double getBase() {
		return this.base;
	}
	
	public double doubleValue() {
		return this.base * Math.pow(2, this.scale);
	}
	
	public SmallDecimal multiply(SmallDecimal factor) {
		double b = this.base * factor.base;
		long s = this.scale + factor.scale;
		if (Math.abs(b)>2) {
			b = b/2; scale++;
		}	
		return new SmallDecimal(b,s);
	}
	
	public SmallDecimal divide(SmallDecimal divisor) {
		if (0==divisor.base) {
			throw new RuntimeException("Division by zero.");
		} else {
			double b = this.base / divisor.base;
			long s = this.scale - divisor.scale;
			if (this.base<divisor.base) {
				b = b/2; scale++;
			}
			return new SmallDecimal(b,s);
		}	
	}
	
	public SmallDecimal add(SmallDecimal summand) {
		double b = 0; long s = 0;
		if (this.scale>=summand.scale) {
			s = this.scale;
			b = this.base + summand.base*Math.pow(2, summand.scale-this.scale);
		} else {
			s = summand.scale;
			b = summand.base + this.base*Math.pow(2, this.scale-summand.scale);
		}
		if (Math.abs(b)>2) {
			b=b/2; scale++;
		}
		return new SmallDecimal(b,s);
	}
	
	public SmallDecimal subtract(SmallDecimal summand) {
		double b = 0; long s = 0;
		if (this.scale>=summand.scale) {
			s = this.scale;
			b = this.base - summand.base*Math.pow(2, summand.scale-this.scale);
		} else {
			s = summand.scale;
			b = this.base*Math.pow(2, this.scale-summand.scale) - summand.scale;
		}
		if (Math.abs(b)>2) {
			b=b/2; scale++;
		}
		return new SmallDecimal(b,s);
	}
}
