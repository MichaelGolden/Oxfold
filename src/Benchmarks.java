import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ox.osscb.Program;

public class Benchmarks {
	public static void main(String [] args)
	{
		File dataDir = new File("datasets/");
		File outputDir = new File("output/");
		outputDir.mkdir();
		System.out.println(dataDir.list().length);
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
		int datasetno = 0;
		for(StructureData s : experimentalStructures)
		{
			long startNano = System.nanoTime();
			StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, 3.0);
			long endNano = System.nanoTime();
			double elapsedNano = (endNano - startNano)/1000000000.0;
			predictedStructure.time = elapsedNano;
			predictedStructures.add(predictedStructure);

			datasetno++;
			System.out.println("dataset "+datasetno+"\t"+elapsedNano);
			try {
				Benchmarks.saveBenchmarksCSV(new File("results.csv"), experimentalStructures, predictedStructures);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
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
