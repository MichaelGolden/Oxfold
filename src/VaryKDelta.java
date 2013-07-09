import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.Program;

public class VaryKDelta {
	
	public static ArrayList<String> avgMetrics; 
	
	public static void main(String [] args)
	{
		//String outputDirString = "output/";
		avgMetrics = new ArrayList<String>(); 
		

		if(args.length != 2){
			System.err.println("Default: 0 <= K <= 10, 0.1 <= Delta <= 1, with increments of 0.1");
			runAuto(10, 1);
		}
		else{
			Constants.IterationCutOff = BigDecimal.valueOf(Double.valueOf(args[1])); 
			double delta_o = Constants.IterationCutOff.doubleValue(); 
			double kArg = Double.valueOf(args[0]);
			runSpecific(kArg, delta_o);
		}
		
		
	}
	
	private static void runSpecific(double KValue, double dValue){
		String outputDirString = "output/";
		Constants.IterationCutOff = BigDecimal.valueOf(Double.valueOf(dValue)); 
		double delta = Constants.IterationCutOff.doubleValue(); 
		if(KValue == 0){
			System.err.println("Cannot use K = 0, use K > 0 instead");
			KValue = 0.1; 
		}
		File dataDir = new File("datasets/");
		File outputDir = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue +"/");
		outputDir.mkdir();
		 
		
		ArrayList<StructureData> experimentalStructures = new ArrayList<>();
		for(File experimentalFile : dataDir.listFiles()) //to do: add 40 file max 
		{
			if(experimentalFile.getName().equals("TestRNAData49.dat")){
				try {
					experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			
		}
		
		Collections.sort(experimentalStructures);
		
		ArrayList<StructureData> predictedStructures = new ArrayList<>();
		
		int maxIter = 1; 
		for(int j  = 0; j < maxIter; j++) //go through select # of datasets
		//for(StructureData s : experimentalStructures) //go through all
		{
			StructureData s = experimentalStructures.get(0);
			//StructureData s = experimentalStructures.get(j);
			predictedStructures.add(VaryKDelta.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, KValue));
			
			try {
				//Don't want to save full CSV for this time, leave in possibility
				
				
				//VaryKDelta.saveBenchmarksCSV(new File("output/" + "K_" + KValue + "/results_K_" + KValue+ ".csv"), experimentalStructures, predictedStructures);
				VaryKDelta.saveBenchmarkAvgCSV( experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
				
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} //end dataset iterations
		
		try{
			//set up for avg
			double sensitivity = 0; 
			double ppv = 0;
			double fscore = 0; 
			double mountainSim = 0; 
			File outputCsvDir = new File(outputDirString + "csv/" );
			outputCsvDir.mkdirs();
			FileWriter writerK = new FileWriter(outputDirString + "csv/results_K_" + KValue + "_Delta_" + dValue + ".csv"); 
			for(String str: avgMetrics) {
			  writerK.write(str);
			  String[] line = str.split(",");
			  sensitivity += Double.valueOf(line[2]);
			  ppv += Double.valueOf(line[3]);
			  fscore += Double.valueOf(line[4]);
			  mountainSim += Double.valueOf(line[5]);
			}
			writerK.write("Average,K_" + KValue + "_Delta_" + dValue +","+ sensitivity/avgMetrics.size()
					+","+ ppv/avgMetrics.size()+","+ fscore/avgMetrics.size()+","+ mountainSim/avgMetrics.size());
			writerK.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
					
			
			
	}
	
	private static  void runAuto(double k, double delta){
		String outputDirString = "output/";
		//variables for K or delta and their ranges. Then move below out to new method. 
				for(double KValue = 0; KValue <= k; KValue += 0.1){ //change this later
					for(double dValue = 0; dValue <= delta; dValue += 0.1){
						Constants.IterationCutOff = BigDecimal.valueOf(Double.valueOf(dValue)); 
						delta = Constants.IterationCutOff.doubleValue(); 
						if(KValue == 0){
							KValue = 0.1; 
						}
						File dataDir = new File("datasets/");
						File outputDir = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue +"/");
						outputDir.mkdir();
						
						ArrayList<StructureData> experimentalStructures = new ArrayList<>();
						for(File experimentalFile : dataDir.listFiles()) //to do: add 40 file max 
						{
							if(experimentalFile.getName().equals("TestRNAData49.dat")){
								try {
									experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}

							
						}
						
						Collections.sort(experimentalStructures);
						

						ArrayList<StructureData> predictedStructures = new ArrayList<>();
						
						int maxIter = 1; 
						for(int j  = 0; j < maxIter; j++) //go through select # of datasets
						//for(StructureData s : experimentalStructures) //go through all
						{
							StructureData s = experimentalStructures.get(0);
							//StructureData s = experimentalStructures.get(j);
							predictedStructures.add(VaryKDelta.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, KValue));
							
							try {
								//Don't want to save full CSV for this time, leave in possibility
								
								
								//VaryKDelta.saveBenchmarksCSV(new File("output/" + "K_" + KValue + "/results_K_" + KValue+ ".csv"), experimentalStructures, predictedStructures);
								VaryKDelta.saveBenchmarkAvgCSV( experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
								
								
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						} //end iteration loop
					}
					
					
				} //end K value
				
				try{
					FileWriter writerK = new FileWriter("output/results_K_Delta.csv"); 
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
			System.out.println(KValue +","+ delta + ","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
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
		//System.out.println("Save benchmarks to "+csvFile.getAbsolutePath());
	}
}
