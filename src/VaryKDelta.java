/*import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.Program;

public class VaryKDelta {
	
	public static ArrayList<String> avgMetrics; 
	public static void main(String [] args)
	{
		String outputDirString = "output/";
		avgMetrics = new ArrayList<String>(); 
		
		double KValue; 
		double delta = Constants.IterationCutOff.doubleValue(); 
		//variables for K or delta and their ranges. Then move below out to new method. 
		for(KValue = 0; KValue <=10; KValue += 0.5){ //change this later
			if(KValue == 0){
				KValue = 0.1; 
			}
			File dataDir = new File("datasets/");
			File outputDir = new File(outputDirString + "K_" + KValue + "/");
			outputDir.mkdir();
			
			ArrayList<StructureData> experimentalStructures = new ArrayList<>();
			for(File experimentalFile : dataDir.listFiles()) //to do: add 40 file max 
			{
				try {
					experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			Collections.sort(experimentalStructures);
			

			ArrayList<StructureData> predictedStructures = new ArrayList<>();
			
			int maxIter = 1; 
			for(int j  = 0; j < maxIter; j++) //go through select # of datasets
			//for(StructureData s : experimentalStructures) //go through all
			{
				StructureData s = experimentalStructures.get(j);
				predictedStructures.add(VaryKDelta.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, KValue));
				
				try {
					//Don't want to save full CSV for this time, leave in possibility
					
					
					VaryKDelta.saveBenchmarksCSV(new File("output/" + "K_" + KValue + "/results_K_" + KValue+ ".csv"), experimentalStructures, predictedStructures);
					//VaryKDelta.saveBenchmarkAvgCSV( experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
					
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			} //end iteration loop
			
		} //end K value
		
		try{
			FileWriter writerK = new FileWriter("output/results_K.csv"); 
			for(String str: avgMetrics) {
			  writerK.write(str);
			}
			writerK.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");

		int [] pairedSites = null;
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
		
			
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File("doc/ppfold.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/ppfold.parameters").getAbsolutePath(),
					"--tree="+outNewick.getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
		}
		else
		{
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
					//"--grammar="+new File("doc/kh.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
		}	

		StructureData structureData = new StructureData(pairedSites);
		return structureData;		
	}
	
	
	public static void saveBenchmarkAvgCSV( List<StructureData> experimentalStructures, List<StructureData> predictedStructures, double KValue, double delta, boolean isLast) throws IOException
	{
		//MODIFY TO SAVE AVERAGE OVER ALL ON SAME RUN
		//BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile));
		
		//writer.write("Name,Sensitivity,PPV,FScore,MountainSim\n");
		double senstivity = 0; 
		double ppv = 0;
		double fscore = 0; 
		double mountainSim = 0; 
		double predictedStructuresSize = predictedStructures.size();
		System.out.println("size" + predictedStructuresSize);
		for(int i = 0 ; i < predictedStructuresSize ; i++)
		{
			StructureData expStructure = experimentalStructures.get(i);
			StructureData predictedStructure = predictedStructures.get(i);
			
			senstivity = BasePairMetrics.calculateSensitivity(expStructure.pairedSites, predictedStructure.pairedSites);
			ppv = BasePairMetrics.calculatePPV(expStructure.pairedSites, predictedStructure.pairedSites);
			fscore = BasePairMetrics.calculateFScore(expStructure.pairedSites, predictedStructure.pairedSites);
			mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(expStructure.pairedSites, predictedStructure.pairedSites);
			
			avgMetrics.add(KValue +","+ delta + ","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
		}
		//fix this
		//if(isLast){
		//	avgMetrics.add(KValue +","+ delta + ","+senstivity/predictedStructuresSize+","+ppv/predictedStructuresSize+
		//		","+fscore/predictedStructuresSize+","+mountainSim/predictedStructuresSize+"\n");

		//}
		//avgMetrics.add(KValue +","+ delta + ","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
		
		//writer.write(KValue +","+ delta + ","+senstivity/predictedStructuresSize+","+ppv/predictedStructuresSize+
		//		","+fscore/predictedStructuresSize+","+mountainSim/predictedStructuresSize+"\n");
		//writer.close();
		
	}
	
	public static void saveBenchmarksCSV(File csvFile, List<StructureData> experimentalStructures, List<StructureData> predictedStructures) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile));
		
		writer.write("Name,Length,Sensitivity,PPV,FScore,MountainSim\n");
		
		for(int i = 0 ; i < predictedStructures.size() ; i++)
		{
			StructureData expStructure = experimentalStructures.get(i);
			StructureData predictedStructure = predictedStructures.get(i);
			
			double senstivity = BasePairMetrics.calculateSensitivity(expStructure.pairedSites, predictedStructure.pairedSites);
			double ppv = BasePairMetrics.calculatePPV(expStructure.pairedSites, predictedStructure.pairedSites);
			double fscore = BasePairMetrics.calculateFScore(expStructure.pairedSites, predictedStructure.pairedSites);
			double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(expStructure.pairedSites, predictedStructure.pairedSites);
			
			writer.write(expStructure.file.getName()+","+expStructure.pairedSites.length+","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
		}
		
		writer.close();
		System.out.println("Save benchmarks to "+csvFile.getAbsolutePath());
	}
}
*/

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.Program;

public class VaryKDelta {
	public static ArrayList<String> avgMetrics; 
	

	public static void main(String [] args)
	{
		String outputDirString = "output/";
		avgMetrics = new ArrayList<String>(); 
		
		double KValue; 
		double delta = Constants.IterationCutOff.doubleValue(); 
		//variables for K or delta and their ranges. Then move below out to new method. 
		for(KValue = 0; KValue <=10; KValue += 0.5){ //change this later
			if(KValue == 0){
				KValue = 0.1; 
			}
			File dataDir = new File("datasets/");
			File outputDir = new File(outputDirString + "K_" + KValue + "/");
			outputDir.mkdir();
			ArrayList<StructureData> experimentalStructures = new ArrayList<>();
			for(File experimentalFile : dataDir.listFiles())
			{
				try {
					experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			Collections.sort(experimentalStructures);
		
		
		ArrayList<StructureData> predictedStructures = new ArrayList<>();
		int maxIter = 5; 
		for(int j  = 0; j < maxIter; j++) //go through select # of datasets
		//for(StructureData s : experimentalStructures) //go through all
		{
			StructureData s = experimentalStructures.get(j);
		
			predictedStructures.add(Benchmarks.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, 3.0));
			
			try {
				Benchmarks.saveBenchmarksCSV(new File("results.csv"), experimentalStructures, predictedStructures);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		}
		
		
		//System.out.println(dataDir.listFiles().length);
	}
	
	public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");

		int [] pairedSites = null;
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
		
			
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File("doc/ppfold.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/ppfold.parameters").getAbsolutePath(),
					"--tree="+outNewick.getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
		}
		else
		{
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
					//"--grammar="+new File("doc/kh.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
		}	

		StructureData structureData = new StructureData(pairedSites);
		return structureData;		
	}
	
	public static void saveBenchmarksCSV(File csvFile, List<StructureData> experimentalStructures, List<StructureData> predictedStructures) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile));
		
		writer.write("Name,Length,Sensitivity,PPV,FScore,MountainSim\n");
		
		for(int i = 0 ; i < predictedStructures.size() ; i++)
		{
			StructureData expStructure = experimentalStructures.get(i);
			StructureData predictedStructure = predictedStructures.get(i);
			
			double senstivity = BasePairMetrics.calculateSensitivity(expStructure.pairedSites, predictedStructure.pairedSites);
			double ppv = BasePairMetrics.calculatePPV(expStructure.pairedSites, predictedStructure.pairedSites);
			double fscore = BasePairMetrics.calculateFScore(expStructure.pairedSites, predictedStructure.pairedSites);
			double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(expStructure.pairedSites, predictedStructure.pairedSites);
			
			writer.write(expStructure.file.getName()+","+expStructure.pairedSites.length+","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
		}
		
		writer.close();
		System.out.println("Save benchmarks to "+csvFile.getAbsolutePath());
	}
}
