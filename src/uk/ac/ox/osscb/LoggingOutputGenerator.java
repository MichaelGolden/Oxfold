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

	public String dumpStructure(int[] structure) {
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
}
