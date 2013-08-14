package uk.ac.ox.osscb;
import java.util.InputMismatchException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class LoggingOutputGenerator implements OutputGenerator{
	final Logger log = LoggerFactory.getLogger(LoggingOutputGenerator.class);

	/**
	 * Output dot-bracket structure into a log file.
	 */
	public void generate(int[] structure) {
		
		if(log.isInfoEnabled()){		
			log.info(dumpStructure(structure));
		}
	}
	
	/**
	 * Rarer, more important messages.
	 * Invokes {@link ProgramOutput#outMsg2(String, Object...)}
	 *  
	 * @param structure
	 */
	public void generateFinal(int[] structure) {
		ProgramOutput.outMsg2(String.format("Final structure is:%n%s", dumpStructure(structure)));
	}

	public static String dumpStructure(int[] structure) {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<structure.length; i++)
		{
			char currentIdx = '+';
			if(structure[i] == -1)
				currentIdx = '.';
			else if(structure[i] > i)
				currentIdx = '(';
			else if(structure[i] < i)
				currentIdx = ')';
			else
				throw new InputMismatchException("Input was not an int array representing a structure");
			sb.append(currentIdx);
		}	
		String string = sb.toString();
		return string;
	}
	


	public static String dumpStructure(Structure structure) {
		StringBuilder sb = new StringBuilder();
		int[] pairings = structure.getPairings();
		boolean[] keepPairs = structure.getKeepPairs();
		for(int i=0; i<pairings.length; i++)
		{
			char currentIdx = '+';
			if(pairings[i] == -1 )//|| keepPairs[i] == false)
				currentIdx = '.';
			else if(pairings[i] > i)// && keepPairs[i])
				currentIdx = '(';
			else if(pairings[i] < i )//&& keepPairs[i])
				currentIdx = ')';
			else
				throw new InputMismatchException("Input was not an int array representing a structure");
			sb.append(currentIdx);
		}	
		String string = sb.toString();
		return string;
	}

	//@Override
	public void generate(Structure structure) {
		// TODO Auto-generated method stub
		if(log.isInfoEnabled()){		
			log.info(dumpStructure(structure));
		}
	}

	//@Override
	public void generateFinal(Structure structure) {
		// TODO Auto-generated method stub
		ProgramOutput.outMsg2(String.format("Final structure is:%n%s", dumpStructure(structure)));
	}
}
