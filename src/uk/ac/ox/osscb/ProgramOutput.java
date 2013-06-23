package uk.ac.ox.osscb;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * System-wide logging, instead of {@link System#out}.
 * The reason is this class uses logback loggers
 * which can be configured w/o program recompilation.
 * 
 * The class has 2 logging methods:
 * <ul>
 * <li>{@link #outMsg(String)} - for 'regular' program-wide messages,
 * can be invoked frequently</li>
 * <li>{@link #outMsg2(String)} - for much less frequent program-wide messages,
 * like 'program started' and 'program finished' messages</li>
 * </ul>
 * The former is typically configured to output into the main program log file
 * while the latter - both the file and the console. Configuration can be changed
 * easily.  
 * 
 * Please be aware that this class is not intended to be used for logging
 * purposes, but for important messages (like program output, etc). 
 * 
 * @author Vladimir
 *
 */
public class ProgramOutput {
	
	private static final Logger progOut = LoggerFactory.getLogger(Constants.ProgramOutputLogName);

	private static final Logger progOut2 = LoggerFactory.getLogger(Constants.ProgramOutput2LogName);

	/**
	 * Output normal, system-wide messages. Normally they don't go to a console.
	 * Exact behaviour depends on logback settings.
	 * 
	 * @param msg
	 */
	public static void outMsg(String msg, Object... args){
		if(progOut.isInfoEnabled()){

			if(null != args && args.length > 0){
				msg = String.format(msg, args);
			}
			progOut.info(msg);
		}
	}
	
	/**
	 * Output rarer (and more important) messages which are supposed to get to console.
	 * Exact behaviour depends on logback settings.
	 * 
	 * @param msg
	 */
	public static void outMsg2(String msg, Object... args){
		if(progOut2.isInfoEnabled()){
			
			if(null != args && args.length > 0){
				msg = String.format(msg, args);
			}
			progOut2.info(msg);
		}
	}
}
