package uk.ac.ox.osscb.analysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ox.osscb.CoFoldAnalogue;
import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.HaasMEAFold;
import uk.ac.ox.osscb.Program;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.visualisation.DataVisualiser;
import uk.ac.ox.osscb.visualisation.SVG;

public class BenchmarksAll {

	public enum Method {
		Oxfold, HaasMEA;
		
		public String toString()
		{
			switch(this)
			{
			case Oxfold:
				return "Oxfold";
			case HaasMEA:
				return "HaasMEA";
			}
			
			return "";
		}
	};
	
	public static void main(String [] args)
	{		
		File dataDir = new File("datasets_5seq/"); // name of the input directory (containing .dat files)
		boolean first_sequence_only = false; // if true, does a single sequence prediction (from the first sequence in the alignment)
		boolean run_visualisations = false; // if true save the structure visualisations (much slower)
		int run_shortest_n_datasets = 45; // will run the shortest n datasets
		
		//Method benchmarkMethod = Method.Oxfold; // change this to modify which method is benchmarked
		Method benchmarkMethod = Method.HaasMEA; // change this to modify which method is benchmarked		
		double weightHaasMEA  = 0.2; // helix formation delta criterion for HaasMEA
		double tauHaasMea = Double.POSITIVE_INFINITY; // if +infinity, doesn't use shortcut weighting function
		
		File outputDir = new File(dataDir+"_output");
		outputDir.mkdir();
		File resultsFile = new File(outputDir.toString()+File.separatorChar+"benchmarks.csv");
	
		boolean startFile = true;
		
		System.out.println(dataDir.list().length);
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for(File experimentalFile : dataDir.listFiles())
		{
			try {
				StructureData data = StructureData.readExperimentalStructureData(experimentalFile);
				
				if(first_sequence_only)
				{
					while(data.sequences.size() > 1)
					{
						data.sequences.remove(data.sequences.size()-1);
						data.sequenceNames.remove(data.sequenceNames.size()-1);
					}
				}
				experimentalStructures.add(data);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		Collections.sort(experimentalStructures);
		

		int m = 0;
		for(StructureData e : experimentalStructures)
		{
				System.out.println(m+"\t"+e.file+"\t"+e.sequences.get(0).length());
				m++;
		}
		
		ArrayList<StructureData> predictedStructures = new ArrayList<StructureData>();
		int datasetno = 0;
		long startNanoTime = System.nanoTime();

		int start = 0;
		int end = run_shortest_n_datasets;
		for(int i = 0 ; i < Math.min(end, experimentalStructures.size()) ; i++)
		{
			StructureData s = experimentalStructures.get(i);
			if(datasetno < start)
			{
				datasetno++;
				predictedStructures.add(null);
				continue;
			}
			long startNano = System.nanoTime();
			
			try {
				IO.writeLine(new File(outputDir+"/"+ s.file.getName()+".experimental"), RNAFoldingTools.getDotBracketStringFromPairedSites(s.pairedSites), true, false);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			StructureData predictedStructure = null;
			switch(benchmarkMethod)
			{
				case Oxfold:
					double weightOxfold = 0.5; // this parameter no longer makes a difference, because the code will decide
					predictedStructure  = BenchmarksAll.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,weightOxfold, false);
					break;
				case HaasMEA:
					predictedStructure = BenchmarksAll.foldHaasMEA(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, weightHaasMEA,tauHaasMea, false);
					break;
			}
			StructureData predictedStructure2 = BenchmarksAll.foldRNAalifold(outputDir, s.file.getName()+"_rnaalifold", s.sequences, s.sequenceNames);
			StructureData predictedStructure3 = PPfold.fold(outputDir, s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true);
			
			System.out.println("Entropy: "+predictedStructure3.entropy+"\t"+predictedStructure3.normalisedEntropy);

			long endNanoTime = System.nanoTime();
			double elapsed = (endNanoTime-startNanoTime)/1e9;
			System.out.println("Elapsed\t"+i+"\t"+elapsed+"s");

			StructureData s1 = predictedStructure;
			StructureData s2 = predictedStructure2;
			
			if(run_visualisations) // if true, create the visualisations
			{
				DataVisualiser visualiser = new DataVisualiser();
				
				s2.sequences = s.sequences;
				s2.sequenceNames = s.sequenceNames;
				s1.sequences = s2.sequences;
				s1.sequenceNames = s2.sequenceNames;
				s.title = "Experimental";
				s1.title = benchmarkMethod.toString();
				s2.title = "PPfold";
	
				SVG full = visualiser.drawComparisonPredicted(s1, s2, s);

				try {
					full.savePNG(new File(outputDir.getAbsolutePath()+File.separatorChar+s.file.getName()+".svg"), new File(outputDir.getAbsolutePath()+File.separatorChar+s.file.getName()+".png"));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
					
			long endNano = System.nanoTime();
			double elapsedNano = (endNano - startNano)/1000000000.0;
			predictedStructure.time = elapsedNano;
			predictedStructures.add(predictedStructure);
		
			datasetno++;
			System.out.println("dataset "+datasetno+"\t"+elapsedNano);
			try {
				if(startFile)
				{
					IO.clearFile(resultsFile);
					IO.writeLine(resultsFile, benchmarkMethod.toString()+",,,,,,RNAalifold,,,,,,PPfold,,,,,", true, true);
					IO.writeLine(resultsFile, "Name,Length,Sensitivity,PPV,FScore,MountainSim,Name,Length,Sensitivity,PPV,FScore,MountainSim,Name,Length,Sensitivity,PPV,FScore,MountainSim,N", true, true);
					startFile = false;
				}
				IO.writeLine(resultsFile, getBenchmarkString(s, s1)+","+getBenchmarkString(s, s2)+","+getBenchmarkString(s, predictedStructure3)+","+s.sequences.size(), true, true);
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight, boolean reverseGrammar)
	{
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+"_"+sequences.size()+".nwk");
		File basePairProbFile = null;
		
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
			String grammarPath = reverseGrammar ? "doc/ppfold_reverse.grammar" : "doc/ppfold.grammar";
			
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File(grammarPath).getAbsolutePath(),
					"--grammar-params="+new File("doc/ppfold.parameters").getAbsolutePath(),
					"--tree="+outNewick.getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
			
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".evol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
		}
		else
		{
			String grammarPath = reverseGrammar ? "doc/kh_reverse.grammar" : "doc/kh.grammar";
			
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					//"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
					"--grammar="+new File(grammarPath).getAbsolutePath(),
					"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath(),
					"--weight="+weight};			
	
			new Program().run(argsArray);
		
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".noevol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
		}	

		StructureData structureData = new StructureData(pairedSites);
		if(basePairProbFile.exists())
		{
			structureData.basePairProbFile =  basePairProbFile;
		}
		return structureData;		
	}
	
	public static StructureData foldCofold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double alpha, double tau, boolean reverseGrammar)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+"_"+sequences.size()+".nwk");
		
		File basePairProbFile = null;

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
		
			
			/*String [] argsArray = {fastaFile.getAbsolutePath(), 
					"--grammar="+new File("doc/ppfold.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/ppfold.parameters").getAbsolutePath(),
					"--tree="+outNewick.getAbsolutePath(),
					"--weight="+weight};	*/		
	
			//new Program().run(argsArray);
			CoFoldAnalogue cofold = new CoFoldAnalogue();
			cofold.foldEvolutionary(fastaFile.getAbsolutePath(), new File(reverseGrammar ? "doc/ppfold_reverse.grammar" : "doc/ppfold.grammar").getAbsolutePath(), new File("doc/ppfold.parameters").getAbsolutePath(), outNewick.getAbsolutePath(), alpha, tau);
			
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".evol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
		}
		else
		{
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					//"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
					"--grammar="+new File(reverseGrammar ? "doc/kh_reverse.grammar" : "doc/kh.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath()};
	
			//CoFoldAnalogue cofold = new CoFoldAnalogue();
			//cofold.foldEvolutionary(fastaFile.getAbsolutePath(), new File("doc/kh.grammar").getAbsolutePath(), new File("doc/kh.parameters").getAbsolutePath(), outNewick.getAbsolutePath(), 0.0, 0.0);
			
			//new Program().run(argsArray);
			System.err.println("Not yet implemented");
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".noevol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
		}	

		StructureData structureData = new StructureData(pairedSites);
		if(basePairProbFile.exists())
		{
			structureData.basePairProbFile =  basePairProbFile;
		}
		return structureData;		
	}
	
	public static StructureData foldHaasMEA(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight, double tau, boolean reverseGrammar)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+"_"+sequences.size()+".nwk");
		File basePairProbFile = null;

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
	
			//new Program().run(argsArray);
			HaasMEAFold haasMEA = new HaasMEAFold();
			haasMEA.foldEvolutionary(fastaFile.getAbsolutePath(), new File(reverseGrammar ? "doc/ppfold_reverse.grammar" : "doc/ppfold.grammar").getAbsolutePath(), new File("doc/ppfold.parameters").getAbsolutePath(), outNewick.getAbsolutePath(), weight, tau);
			
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".evol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
		}
		else
		{
			String [] argsArray = {fastaFile.getAbsolutePath(), 
					//"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
					"--grammar="+new File(reverseGrammar ? "doc/kh_reverse.grammar" : "doc/kh.grammar").getAbsolutePath(),
					"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath()};
	
			//CoFoldAnalogue cofold = new CoFoldAnalogue();
			//cofold.foldEvolutionary(fastaFile.getAbsolutePath(), new File("doc/kh.grammar").getAbsolutePath(), new File("doc/kh.parameters").getAbsolutePath(), outNewick.getAbsolutePath(), 0.0, 0.0);
			
			//new Program().run(argsArray);
			System.err.println("Not yet implemented");
			basePairProbFile = new File(fastaFile.getAbsolutePath()+".noevol.bp");
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
		}	

		StructureData structureData = new StructureData(pairedSites);
		if(basePairProbFile.exists())
		{
			structureData.basePairProbFile =  basePairProbFile;
		}
		return structureData;		
	}
	

	
	public static StructureData foldRNAalifold(File dir, String name, ArrayList<String> sequences, ArrayList<String> sequenceNames)
	{
		StructureData data = new StructureData();

		try {
			RNAalifoldResult res = RNAalifold.fold(sequences, sequenceNames, "", true, false);
			data.sequences = sequences;
			data.sequenceNames = sequenceNames;
			data.pairedSites = res.pairedSites;
			File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
			File basePairProbFile = new File(fastaFile.getAbsolutePath()+".bp");			
			if(data.basePairProb != null)
			{
				RNAFoldingTools.saveMatrix(basePairProbFile, res.matrix);	
				data.basePairProbFile = basePairProbFile;
			}
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		return data;
	}
	
	public static String getDistributionString(StructureData data)
	{
		try
		{
			ArrayList<Double> all = new ArrayList<Double>();
			ArrayList<Double> structure = new ArrayList<Double>();
			double [][] bp = StructureData.getBasePairProb(data.basePairProbFile);
			double [][] matrix = PosteriorProbabilitiesCalculator.getDiffs(bp);
			for(int i = 0 ; i < matrix.length ; i++)
			{
				for(int j = i+1 ; j < matrix.length ; j++)
				{
					all.add(matrix[i][j]);
					
				}
			}
			
			for(int i = 0 ; i < data.pairedSites.length ; i++)
			{
				if(data.pairedSites[i] != 0 && i < data.pairedSites[i]-1)
				{
					structure.add(matrix[i][data.pairedSites[i]-1]);
				}
			}
			
			System.out.println("all size "+all.size()+"\t"+bp.length);
			System.out.println("all size "+structure.size());
		
			double [] allDist = {RankingAnalyses.getMin(all), RankingAnalyses.getPercentile(all, 0.25), RankingAnalyses.getMedian(all), RankingAnalyses.getPercentile(all, 0.75), RankingAnalyses.getMax(all)};
			double [] structureDist = {RankingAnalyses.getMin(structure), RankingAnalyses.getPercentile(structure, 0.25), RankingAnalyses.getMedian(structure), RankingAnalyses.getPercentile(structure, 0.75), RankingAnalyses.getMax(structure)};
		
			String ret = "";
			for(int i = 0 ; i < allDist.length ; i++)
			{
				ret += allDist[i]+",";
			}
			for(int i = 0 ; i < structureDist.length ; i++)
			{
				ret += structureDist[i];
				if(i != structureDist.length - 1)
				{
					ret += ",";
				}
			}
			
			return ret;
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		
		return null;
	}
	
	
	public static String getBenchmarkString(StructureData expStructure, StructureData predictedStructure)
	{
		double senstivity = BasePairMetrics.calculateSensitivity(expStructure.pairedSites, predictedStructure.pairedSites);
		double ppv = BasePairMetrics.calculatePPV(expStructure.pairedSites, predictedStructure.pairedSites);
		double fscore = BasePairMetrics.calculateFScore(expStructure.pairedSites, predictedStructure.pairedSites);
		double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(expStructure.pairedSites, predictedStructure.pairedSites);
		
		return expStructure.file.getName()+","+expStructure.pairedSites.length+","+senstivity+","+ppv+","+fscore+","+mountainSim;
	}
	
	public static void saveBenchmarksCSV(File csvFile, List<StructureData> experimentalStructures, List<StructureData> predictedStructures) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile));
		
		writer.write("Name,Length,Sensitivity,PPV,FScore,MountainSim\n");
		
		for(int i = 0 ; i < predictedStructures.size() ; i++)
		{
			StructureData expStructure = experimentalStructures.get(i);
			StructureData predictedStructure = predictedStructures.get(i);
			
			if(predictedStructure != null)
			{
				double senstivity = BasePairMetrics.calculateSensitivity(expStructure.pairedSites, predictedStructure.pairedSites);
				double ppv = BasePairMetrics.calculatePPV(expStructure.pairedSites, predictedStructure.pairedSites);
				double fscore = BasePairMetrics.calculateFScore(expStructure.pairedSites, predictedStructure.pairedSites);
				double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(expStructure.pairedSites, predictedStructure.pairedSites);
				
				writer.write(expStructure.file.getName()+","+expStructure.pairedSites.length+","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
			}
		}
		
		writer.close();
		System.out.println("Save benchmarks to "+csvFile.getAbsolutePath());
	}
}
