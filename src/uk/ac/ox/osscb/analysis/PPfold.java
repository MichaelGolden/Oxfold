package uk.ac.ox.osscb.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class PPfold {
	public static String PPFOLD_EXECUTABLE = "binaries/PPfold3.0.jar";
	
	public static StructureData fold(File dir, String name, ArrayList<String> sequences, ArrayList<String> sequenceNames, boolean runEvolutionary)
	{
		return fold(dir,name,sequences,sequenceNames,runEvolutionary,false);
	}

	public static StructureData fold(File dir, String name, ArrayList<String> sequences, ArrayList<String> sequenceNames, boolean runEvolutionary, boolean useCache)
	{	
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");
		File structureFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".ct");
		File entropyFile =  new File(dir.getAbsolutePath()+File.separatorChar+name+".entropy");
		
		
		if(structureFile.exists())
		{
			StructureData d = new StructureData(RNAFoldingTools.getPairedSitesFromCtFile(structureFile));
			d.basePairProbFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".bp");
			System.out.println("Using cached");
			try
			{
				BufferedReader buffer = new BufferedReader(new FileReader(entropyFile));
				d.entropy = Double.parseDouble(buffer.readLine());
				d.normalisedEntropy = Double.parseDouble(buffer.readLine());
				buffer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}
			
			return d;
		}

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


        double entropy = 0;
		try {
			Process process = Runtime.getRuntime().exec(cmd, null, dir);

			//FastTree.nullOutput(process.getInputStream());
			//FastTree.nullOutput(process.getErrorStream());
			
			BufferedReader r = new BufferedReader(new InputStreamReader(process.getInputStream()));
	          String textline = null;
	          while((textline = r.readLine()) != null)
	          {
	        	 if(textline.matches("ENTROPY:.+=\\s.*"))
	        	 {
	        		 entropy = Double.parseDouble(textline.replaceAll("ENTROPY:.+=\\s", "").trim());	        		
	        	 }
	        		 
	          }

			int exitCode = process.waitFor();
			System.out.println(exitCode);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		StructureData d = new StructureData(RNAFoldingTools.getPairedSitesFromCtFile(structureFile));
		d.basePairProbFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".bp");
		d.entropy = entropy;
		double maxentropy_nr = 0.142- 1.5*Math.log(d.pairedSites.length)/Math.log(2) + 1.388*d.pairedSites.length;
		//double maxentropy = new Double(df2.format(maxentropy_nr)).doubleValue();
		d.normalisedEntropy = (entropy/maxentropy_nr);
		try {
			IO.writeLine(entropyFile, d.entropy+"\n"+d.normalisedEntropy, true, false);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return d;
	}
}
