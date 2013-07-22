package uk.ac.ox.osscb.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class PPfold {
	public static String PPFOLD_EXECUTABLE = "binaries/PPfold3.0.jar";

	public static StructureData fold(File dir, String name, ArrayList<String> sequences, ArrayList<String> sequenceNames, boolean runEvolutionary)
	{	
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");

		if(runEvolutionary)
		{
			if(!outNewick.exists())
			{
				try {
					FastTree.nucleotideGTRFastTree(fastaFile, outNewick);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

		//String cmd = "java -jar " + new File(PPFOLD_EXECUTABLE).getAbsolutePath() + " \"" + fastaFile.getAbsolutePath() + "\" --onlyCT ";
		String cmd = "java -jar " + new File(PPFOLD_EXECUTABLE).getAbsolutePath() + " \"" + fastaFile.getAbsolutePath() + "\" --exports ";
		if(runEvolutionary)
		{
			cmd += "-t \""+outNewick.getAbsolutePath()+"\"";
		}  

		try {
			Process process = Runtime.getRuntime().exec(cmd, null, dir);
			FastTree.nullOutput(process.getInputStream());
			FastTree.nullOutput(process.getErrorStream());
			
			/*BufferedReader r = new BufferedReader(new InputStreamReader(process.getInputStream()));
	          String textline = null;
	          while((textline = r.readLine()) != null)
	          {
	        	//  System.out.println("error: "+textline);
	          }*/

			int exitCode = process.waitFor();
			System.out.println(exitCode);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		StructureData d = new StructureData(RNAFoldingTools.getPairedSitesFromCtFile(new File(dir.getAbsolutePath()+File.separatorChar+name+".ct")));
		d.basePairProbFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".bp");

		return d;
	}
}
