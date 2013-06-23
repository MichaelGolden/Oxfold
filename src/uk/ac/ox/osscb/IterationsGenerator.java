package uk.ac.ox.osscb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class IterationsGenerator {
	private final Logger log = LoggerFactory.getLogger(IterationsGenerator.class);
	private LoggingOutputGenerator loggingOutputGen = new LoggingOutputGenerator();
	private long lastIterTime = -1;
	private long nonLogIters = 0;
	private long totalIters = 0;
	
	private final int sec = 1;
	private final int iterations = 1;
	
	public void generate(){
		generate(null);
	}

	
	public void generate(int[] struct){
		
		this.nonLogIters++;
		this.totalIters++;
		
		if(this.nonLogIters > iterations){
			doOut(struct);
			return;
		}
		
		long currentTimeMillis = System.currentTimeMillis();
		double secPassed = (currentTimeMillis - lastIterTime) / 1000.0;
		
		if(secPassed > sec){
			doOut(struct);
			return;
		}
	}
	
	private void doOut(int[] struct){
		if(log.isInfoEnabled()){
			String iterations = String.format("Iterations: %s", this.totalIters);
			log.info(iterations);
		
			if(null != struct)
				loggingOutputGen.generate(struct);	
		}
		/*
		System.out.println(String.format("Iterations: %s", this.totalIters));
		
		if(null != struct){
			System.out.println(String.format("\tstructure: %s", Arrays.toString(struct)));			
		} */

		this.lastIterTime = System.currentTimeMillis();
		this.nonLogIters = 0;
		
	}
}
