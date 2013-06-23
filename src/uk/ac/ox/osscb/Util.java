package uk.ac.ox.osscb;

import java.io.File;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Date;
import java.util.Locale;

import org.joda.time.Period;
import org.joda.time.format.PeriodFormatter;
import org.joda.time.format.PeriodFormatterBuilder;

public class Util {
	
	/**
	 * Use %n instead of getting system property because the latter may
	 * potentially throw security exceptions as per code signature and it's easier
	 * to use String. Although it's highly probable that %n boils down to getting
	 * system property value.
	 */
	private static final String NEWLINE = String.format("%n");

	public static BigDecimal[][] makeSquareZerosMatrix(int dim){
		
		BigDecimal[][] bigDecAr = new BigDecimal[dim][dim];
		
		for(int i = 0; i < dim; i++){
			for(int j = 0; j < dim; j++){
				bigDecAr[i][j] = BigDecimal.ZERO;
			}
		}
		return bigDecAr;
	}
	
	public static String dump1DArray(BigDecimal[] numAr){
		return dump1DArray(numAr, 3, 5);
	}
	
	public static String dump1DArray(BigDecimal[] numAr, int precision, int width){
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		String formatStr = createFormatStringG(precision, width);
		for(int i = 0; i < numAr.length; i++){
			if(i > 0)
				sb.append(", ");
			BigDecimal val = numAr[i];
			String valStr;
			if(null != val){
				valStr = String.format(Locale.UK, formatStr, val);
			}else{
				valStr = "<nl>";
			}
			sb.append(valStr);
		}
		sb.append("]");
		
		return sb.toString();
	}
	
	public static String print1DArray(double[] numAr, int precision, int width){
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		String formatStr = createFormatStringG(precision, width);
		for(int i = 0; i < numAr.length; i++){
			if(i > 0)
				sb.append(", ");
			String valStr = String.format(Locale.UK, formatStr, numAr[i]);
			sb.append(valStr);
		}
		sb.append("]");
		
		return sb.toString();
	}	
	
	private static String createFormatStringG(int precision, int width) {
		return createFormatString(precision, width, 'g');
	}
	
	private static String createFormatStringF(int precision, int width) {
		return createFormatString(precision, width, 'f');
	}

	private static String createFormatString(int precision, int width, char convertion) {
		return String.format(Locale.UK, "%%%d.%d%c", width, precision, convertion);
	}
	
	public  static String print2DArray(BigDecimal[][] numAr) {
		return print2DArray(numAr, 5, 2);
	}

	/**
	 * TODO make width & precision be 'active'. Currently they're ignored.
	 * 
	 * @param numAr
	 * @param width
	 * @param precision
	 * @return
	 */
	public  static String print2DArray(BigDecimal[][] numAr, int width, int precision) {
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < numAr.length; i++){
			sb.append("[");
			for(int j = 0; j < numAr[i].length;j++){
				BigDecimal val = numAr[i][j];
				String format;
				if(null != val){
					double doubleValue = val.doubleValue();
					format = String.format("%8.2g; ", doubleValue);
				}else{
					format = "<null>";
				}
				sb.append(format);
			}
			sb.append("]").append(nL());
		}
		
		return sb.toString();
	}
	
	public  static String print2DArray(double[][] numAr) {
		return print2DArray(numAr, 5, 2);
	}

	/**
	 * TODO make width & precision be 'active'. Currently they're ignored.
	 * 
	 * @param numAr
	 * @param width
	 * @param precision
	 * @return
	 */
	public  static String print2DArray(double[][] numAr, int width, int precision) {
		StringBuilder sb = new StringBuilder();
		String formatStr = createFormatStringG(precision, width) + "; ";
		for(int i = 0; i < numAr.length; i++){
			sb.append("[");
			for(int j = 0; j < numAr[i].length;j++){
				double val = numAr[i][j];
				String format = String.format(formatStr, val);					
				sb.append(format);
			}
			sb.append("]").append(nL());
		}
		
		return sb.toString();
	}
	
	
	public  static void print2DArray(StringBuilder sb, BigDecimal[][] numAr) {
		for(int i = 0; i < numAr.length; i++){
			sb.append(Arrays.toString(numAr[i]));
			sb.append(nL());
		}
	}
	
	public  static void print2DArray(StringBuilder sb, boolean[][] numAr) {
		for(int i = 0; i < numAr.length; i++){
			sb.append(Arrays.toString(numAr[i]));
			sb.append(nL());
		}
	}
	
	public static String time(){
		return String.format("%1$tH:%1$tM:%1$tS", new Date());
	}
	
	public static String date(){
		return String.format("%1$tY.%1$tm.%1$td", new Date());
	}
	
	public static String dateTime(){
		return String.format("%s %s", date(), time());
	}
	
	/**
	 * To use output of {@link System#nanoTime()}
	 * 
	 * @param nanoSec
	 * @return
	 */
	public static String spentTimeFromNanoSec(long nanoSec){
		
		double d = (double)(nanoSec) / 1e6;
		long rounded = Math.round(d );
		Period period = new Period( rounded).normalizedStandard();
		
		PeriodFormatter daysHoursMinutes = new PeriodFormatterBuilder()
	    .appendDays()	    
	    .appendSuffix(" d")
	    .printZeroAlways()
	    .appendSeparator(" and ")
	    .appendHours()
	    .appendSuffix(":")
	    .appendMinutes()
	    .appendSuffix(":")
	    .appendSeconds()
	    .appendSuffix(".")
	    .appendMillis()
	    //.appendSuffix(" ms", " ms")
	    .toFormatter();
		
		return daysHoursMinutes.print(period);
	}

	/**
	 * 
	 * @param fileName
	 * @return fileName <b>exactly</b> as it was passed to the function.
	 * In order to be able to chain the function invocation.
	 * 
	 */
	public static String assertCanReadFile(String fileName){
		
		if(null == fileName || fileName.trim().length() < 1){
			throw new IllegalArgumentException("fileName cannot be empty");
		}
		
		File file = new File(fileName);
		if (!file.exists()) {
			throw new IllegalArgumentException(String.format("file: '%s' cannot be found", fileName));
		}
		if (!file.canRead()) {
			throw new IllegalArgumentException(String.format("file: '%s' is not readable (not enough permissions?)", fileName));
		}
		if(file.isDirectory()){
			throw new IllegalArgumentException(String.format("file: '%s' is a directory, text file is expected", fileName));			
		}
		
		return fileName;
	}
	
	public static String nL(){
		return NEWLINE;
	}
}
