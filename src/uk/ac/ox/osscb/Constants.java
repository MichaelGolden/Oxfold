package uk.ac.ox.osscb;

import java.math.BigDecimal;

import org.slf4j.MDC;

public class Constants {
	
	public final static int		UnpairedBaseIdx = -1;
	
	public final static double	DefaultWeightParam = 100;
	
	public final static double	MaxWeightExit = 0.3;
	
	/**
	 * Maximum number of iterations within main workhorse 
	 * {@link KineticFold#fold(String, String, String, double)}
	 * function in order not to get into it forever.
	 * 
	 */
	public final static int		MaxIterations = 1000;
	
	/**
	 * Set it to zero to get more iterations
	 */
	public static  Double IterationCutOffDouble = 0.5;
	public  static BigDecimal IterationCutOff = BigDecimal.valueOf(IterationCutOffDouble);
	
	//public static  Double IterationCutOffDouble;
	//public  static BigDecimal IterationCutOff;
	
	/*public Double getIterationCutOffDouble(){
		IterationCutOffDouble = Delta.IterationCutOffDouble;
		return IterationCutOffDouble.doubleValue();
	}
	*/
	/*public void setIterationCutOffDouble(double d){
		IterationCutOffDouble = d;
		IterationCutOff = BigDecimal.valueOf(IterationCutOffDouble);
	}*/
	
	/**
	 * cutoff for switch from heuristic to evolutionary model
	 */
	public final static double EndNonEvolutionaryFoldDouble = 0.99;
	
	public final static BigDecimal EndNonEvolutionaryFold = BigDecimal.valueOf(EndNonEvolutionaryFoldDouble);
	
	/**
	 * cutoff for switch to full dynamic tests
	 */
	public final static BigDecimal DynamicCutOff = BigDecimal.valueOf(0.99);
	
	public final static int TranscriptionOffset = 4;
	
	public final static int TranscriptionJump = 20;
	
	public final static int TranscriptionMinimumHelixLength = 3;
	
	public final static BigDecimal BreakingCutOff = BigDecimal.valueOf(0.2);
	
	/*public final static BigDecimal priorWeight = BigDecimal.valueOf(3);
	
	public final static BigDecimal priorWeightMinusOne = priorWeight.subtract(BigDecimal.ONE);
	*/
	
	public final static BigDecimal MEApenalty = BigDecimal.valueOf(0.4);
	
	/**
	 * Name of the logger to be used to output main program info. 
	 */
	public final static String ProgramOutputLogName = "OxFold.output";
	
	/**
	 * Name of the logger which is supposed to have only very limited
	 * number of messages which are typically (but not necessarily) 
	 * will get to both console and the main program log file.
	 * 
	 */
	public final static String ProgramOutput2LogName = "OxFold.output2";

	/**
	 * key of logback {@link MDC} which is used to specify 'type' of run
	 * which is combined sequnece number plus a letter from {A, B, C, D}
	 */
	public static final String RUN_SUFFIX = "runSuffix";

	/**
	 * 
	 */
	public static final String MAIN_PROPERTIES_FILE = "oxfold.properties";
}
