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
import uk.ac.ox.osscb.Program;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.visualisation.DataVisualiser;
import uk.ac.ox.osscb.visualisation.SVG;

public class Benchmarks {
	public static void main(String [] args)
	{
		File dataDir = new File("datasets/");
		File outputDir = new File("output13/");
		File resultsFile = new File("results_oxfold_evol_test.csv");

		boolean startFile = true;
		
		outputDir.mkdir();
		System.out.println(dataDir.list().length);
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for(File experimentalFile : dataDir.listFiles())
		{
			try {
				StructureData data = StructureData.readExperimentalStructureData(experimentalFile);
				
				/*while(data.sequences.size() > 1)
				{
					data.sequences.remove(data.sequences.size()-1);
					data.sequenceNames.remove(data.sequenceNames.size()-1);
				}*/
				/*
				for(int i = 0 ; i < 4 ; i++)
				{
					data.sequences.add(data.sequences.get(0));
					data.sequenceNames.add(data.sequenceNames.get(0));
				}*/
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
		//int start = 0;
		//int end = 1;
		int start = 0;
		int end = 45;
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
			
			double oxfoldK = 0.5;
			//StructureData predictedStructure = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0,Double.POSITIVE_INFINITY, false);
			//StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,oxfoldK, false);
			try {
				IO.writeLine(new File(outputDir+"/"+ s.file.getName()+".experimental"), RNAFoldingTools.getDotBracketStringFromPairedSites(s.pairedSites), true, false);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			
			//StructureData predictedStructure = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0.5,500, false);
			//StructureData predictedStructure = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0,250, false);
			//StructureData predictedStructure2 = predictedStructure;
			//System.exit(0);
			//StructureData predictedStructure2 = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0,640, false);
			//StructureData predictedStructure = PPfold.fold(outputDir, s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true);
			//////StructureData predictedStructure2 = PPfold.fold(outputDir, s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true);
			StructureData predictedStructure2 = Benchmarks.foldRNAalifold(outputDir, s.file.getName()+"_rnaalifold", s.sequences, s.sequenceNames);
		//	System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(predictedStructure2.pairedSites));
			//StructureData predictedStructure2 = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			
			
			StructureData predictedStructure3 = PPfold.fold(outputDir, s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true);
			
			System.out.println("Entropy: "+predictedStructure3.entropy+"\t"+predictedStructure3.normalisedEntropy);
			/*
			int [] decodedSites= null;
			try {
				decodedSites=RNAFoldingTools.getPosteriorDecodingConsensusStructure(predictedStructure.getBasePairProb(predictedStructure.basePairProbFile));
				System.err.println(RNAFoldingTools.getDotBracketStringFromPairedSites(decodedSites));
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			StructureData predictedStructure2 = new StructureData(decodedSites);*/
			
			//System.exit(0);
			//System.out.println(s.file.getName());
			//StructureData predictedStructure2 = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0,640, false);
			//System.exit(0);
			//System.exit(0);
			//System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(predictedStructure.pairedSites));
			//System.exit(0);
			//StructureData predictedStructure2 = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, 0,640, false);
			//StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			//StructureData predictedStructure2 = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			
			//StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			//StructureData predictedStructure2 = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			
			
			//StructureData predictedStructure2 = Benchmarks.foldCofold(outputDir, s.file.getName()+"_reverse", s.sequences, s.sequenceNames, true, 0.5,640, false);
		
			long endNanoTime = System.nanoTime();
			double elapsed = (endNanoTime-startNanoTime)/1e9;
			System.out.println("Elapsed\t"+i+"\t"+elapsed+"s");
			
			boolean print = true;
			if(print)
			{
				/*StructureData s1 = predictedStructure;
				s1.sequences = s.sequences;
				s1.sequenceNames = s.sequenceNames;
				StructureData s2 = s;*/
	
				DataVisualiser visualiser = new DataVisualiser();
	
				//StructureData predictedStructure = Benchmarks.foldCofold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, true, 0,640);
				//StructureData predictedStructure2 = Benchmarks.foldCofold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, true, 0.5,640);
				
				StructureData s1 = predictedStructure;
				StructureData s2 = predictedStructure2;
				
				
				s2.sequences = s.sequences;
				s2.sequenceNames = s.sequenceNames;
				s1.sequences = s2.sequences;
				s1.sequenceNames = s2.sequenceNames;
				//s1.title = "Predicted";
				//s2.title = "Experimental";
				s.title = "Experimental";
				//s1.title = "Cofold (null)";
				//s2.title = "Cofold (alpha=0.5 tau=640)";
				s1.title = "Cofold";
				s2.title = "PPfold";
	
				//SVG full = visualiser.drawComparisonPredicted(s1, s2, s);
	
	
	
				//SVG full = visualiser.drawComparisonPredictedExperimental(s1, s2);
				/*try {
					full.savePNG(new File(outputDir.getAbsolutePath()+File.separatorChar+s.file.getName()+".svg"), new File(outputDir.getAbsolutePath()+File.separatorChar+s.file.getName()+".png"));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}*/
				
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
						IO.writeLine(resultsFile, "Oxfold (delta="+Constants.IterationCutOffDouble+" k="+oxfoldK+" gap%="+Constants.gapPercentage+"),,,,,,RNAalifold,,,,,,PPfold,,,,,", true, true);
						IO.writeLine(resultsFile, "Name,Length,Sensitivity,PPV,FScore,MountainSim,Name,Length,Sensitivity,PPV,FScore,MountainSim,Name,Length,Sensitivity,PPV,FScore,MountainSim,N", true, true);
						startFile = false;
					}
					IO.writeLine(resultsFile, getBenchmarkString(s, s1)+","+getBenchmarkString(s, s2)+","+getBenchmarkString(s, predictedStructure3)+","+s.sequences.size(), true, true);
					//Benchmarks.saveBenchmarksCSV(resultsFile, experimentalStructures, predictedStructures);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
		}
		
		
		
		//System.out.println(dataDir.listFiles().length);
	}
	
	public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight, boolean reverseGrammar)
	{
		for(int i = 0 ; i < sequenceNames.size() ; i++)
		{
			sequenceNames.set(i, "seq"+i);
		}
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");
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
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");
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
