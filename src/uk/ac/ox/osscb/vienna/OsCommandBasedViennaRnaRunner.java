package uk.ac.ox.osscb.vienna;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.NoSuchElementException;
import java.util.Scanner;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.ProgramOutput;
import uk.ac.ox.osscb.Util;

/**
 * Invokes Vienna RNAalifold in '-C' mode, like:
 * <i>cat ./OxFoldStruct4.dat | ../../vienna/bin/RNAalifold -p --constraint TestRNAData4.stock.dat</i>
 * Where <i>cat ./OxFoldStruct4.dat</i> is a constraining structure. 
 * 
 * @author Vladimir
 *
 */
public class OsCommandBasedViennaRnaRunner implements ViennaRunner{
	
	private static final Logger log = LoggerFactory.getLogger(OsCommandBasedViennaRnaRunner.class);

	private String path2Vienna;
	private String command;
	private Integer seqLength;
	
	public OsCommandBasedViennaRnaRunner(String path2Vienna, String stockFile, Integer seqLength) {
		super();
		
		Util.assertCanReadFile(stockFile);
		// at least that we can read it.
		// TODO: change for 'x' permission instead
		Util.assertCanReadFile(path2Vienna);
		
		
		this.path2Vienna = path2Vienna;
		this.command = String.format("%s -p --constraint %s", this.path2Vienna, stockFile);
		this.seqLength = seqLength;
	}

	public String getRawOutput(String constrainingStructure/*, boolean useEchoForStdin*/) {
		
		if(StringUtils.isBlank(constrainingStructure)){
			throw new IllegalArgumentException("constrainingStructure must not be empty!");
		}

		String command2;
//		if(useEchoForStdin){
//			// log.warn(String.format)
//			command2 = String.format("echo \"%s\" | %s", constrainingStructure, this.command);
//		}else{
			// we'll use process's output stream to write stdin contents
			command2 = String.format("%s", this.command);
//		}

		log.debug("RNAalifold command:\n\t{}", command2);

		Runtime rt = Runtime.getRuntime();
		try {
			Process process = rt.exec(command2);
			
			//if(!useEchoForStdin){
				// ProgramOutput.outMsg2("\tattempting to write to process output stream: %s", constrainingStructure);
				// do not use echo for stdin, but rather process's own output stream
				OutputStream stdout = null;
				OutputStreamWriter osw = null;
				try{
					stdout = process.getOutputStream();
					osw = new OutputStreamWriter(stdout); 
					osw.write(constrainingStructure);
				}finally{
					// a bit of overkill
					try{
						if(null != osw)
							osw.close();
					}finally{
						if(null != stdout)
							stdout.close();					
					}
				}
			//}

			String programOutput = null;
			Scanner sc = null;
			try{
				sc = new Scanner(process.getInputStream()).useDelimiter("\\A");
				if(sc.hasNext()) {
					programOutput = sc.next();
				}
			}finally{
				if(null != sc)
					sc.close();
			}
			
			return programOutput;
		} catch (IOException ex) {
			throw new RuntimeException(String.format("failed to execute program '%s'. Error: %s",
					this.path2Vienna, ex.getMessage()), ex);
		}
	}

	@Override
	public double[][] getViennaProbabilities(String constrainingStructure) {
		
		getRawOutput(constrainingStructure);
		
		File viennaOutputFile = new File("alidot.ps");
		if(!viennaOutputFile.exists()){
			// provide some help:
			File curDir = new File(".");
			String curDirPath = curDir.getAbsolutePath();
			String[] subFiles = curDir.list();
			StringBuilder sb = new StringBuilder();
			sb.append("Current dir: ").append(curDirPath).append(Util.nL());
			if(null == subFiles){
				sb.append("\tdoes not have any files");				
			}else{
				sb.append("\thas the following files:").append(Util.nL());				
				for(String subFile : subFiles){
					sb.append("\t\t").append(subFile).append(Util.nL());
				}
				sb.append(String.format("\t%d in total", subFiles.length));				
			}
			
			throw new IllegalStateException(String.format("Output file: %s not found. Context:%n%s", 
					viennaOutputFile.getName(), sb));
		}
		
		String next;
		try {
			next = new Scanner(viennaOutputFile).useDelimiter("\\A").next();
		} catch (FileNotFoundException e) {
			// should never get here
			throw new IllegalArgumentException("Couldn't find file: " + viennaOutputFile.getAbsolutePath());
	    } catch (NoSuchElementException e) {
			throw new IllegalArgumentException(String.format("file '%s' is empty", viennaOutputFile.getAbsolutePath()));
	    }

		return new ViennaDotPsFileParser(seqLength).parse(next);
	}
}
