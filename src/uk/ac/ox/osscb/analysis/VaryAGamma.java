package uk.ac.ox.osscb.analysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import uk.ac.ox.osscb.EvolutionaryParameters;
import uk.ac.ox.osscb.EvolutionaryTree;
import uk.ac.ox.osscb.EvolutionaryTreeParser;
import uk.ac.ox.osscb.ParameterParserEvolutionary;
import uk.ac.ox.osscb.Program;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class VaryAGamma {
	
	
	public static void main(String [] args)
	{		
		boolean runEvol = true;

	
		int numdatasets = 25;
		if(args.length > 3)
		{
			numdatasets = Integer.parseInt(args[0]);
			double a = Double.parseDouble(args[1]);
			

			
			
			ArrayList<AGammaJobInput> inputs = new ArrayList<AGammaJobInput>();
			for(int i = 2 ; i < args.length ; i++)
			{
				double gamma = Double.parseDouble(args[i]);
				inputs.add(new AGammaJobInput(a, gamma, runEvol, numdatasets));
			}
			VaryAGamma varyAGamma = new VaryAGamma();
			 try {
				
				List<AGammaJobOutput> outputs = varyAGamma.processInputs(inputs);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		else
		if(args.length != 2 && args.length != 3){
			
			double [] as = new double[39];
			double minA = 1 / (double)(as.length-1);
			for(int i = 0 ; i < as.length ; i++)
			{
				double alpha = -Math.log10(0.15);
				as[i] = Math.pow(10, alpha*(2*(minA*i) - 1));
				System.out.print(as[i]+ " ");
			}
			System.out.println();
			System.exit(0);
			
			double maxgamma = 1;
			double [] gammas = new double[20];
			double mingamma = 0;
			double incgamma = (maxgamma-mingamma)/((double)gammas.length);
			for(int i = 0 ; i < gammas.length ; i++)
			{
				gammas[i] = mingamma + incgamma*i;
				System.out.print(gammas[i]+ " ");
			}
			
			//double [] ks = {0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
			//double [] gammas = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
			

			ArrayList<AGammaJobInput> inputs = new ArrayList<AGammaJobInput>();
			for(double a : as)
			{
				for(double gamma : gammas)
				{
					inputs.add(new AGammaJobInput(a, gamma, runEvol, numdatasets));
				}
			}
			VaryAGamma varyAGamma = new VaryAGamma();
			 try {
				

				List<AGammaJobOutput> outputs = varyAGamma.processInputs(inputs);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			/*
			System.err.println("Default: 0 <= A <= 10, 0.1 <= Gamma <= 1, with increments of 0.1");
			runAuto(10, 1);*/
		}
		else{
			numdatasets = Integer.parseInt(args[2]);
			
			double gamma_o = Double.valueOf(args[1]);
			double aArg = Double.valueOf(args[0]);
			if(args.length > 2)
			{
				runSpecific(aArg, gamma_o, numdatasets, runEvol);
			}
			else
			{
				runSpecific(aArg, gamma_o, numdatasets, runEvol);
			}
		}
		
	}
	
	private static void runSpecific(double AValue, double gValue, int uptodataset, boolean runEvol){
		ArrayList<String> avgMetrics = new ArrayList<String>();
		String outputDirString = "output/";
		runEvol = true; 

		File dataDir = new File("datasets/");
		String evol = runEvol ? "_evol" : "_noevol";
		System.out.println("evol hit 3 = " + evol);
		File outputDir = new File(outputDirString + "A_" + AValue + "_Gamma_" + gValue +evol+"/");
		outputDir.mkdirs();
		 
		
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for(File experimentalFile : dataDir.listFiles())
		{
				try {
					experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}			
		}
		
		Collections.sort(experimentalStructures); // sort datasets in order of increasing length
		

		int maxIter = uptodataset;
		ArrayList<StructureData> sublist =  new ArrayList<StructureData>();
		sublist.addAll((List<StructureData>)experimentalStructures.subList(0, maxIter));
		experimentalStructures = sublist;
		
		ArrayList<StructureData> predictedStructures = new ArrayList<StructureData>();
		 
		int j = 0;
		for(StructureData s : experimentalStructures) 
		{
			//StructureData s = experimentalStructures.get(0);
			//StructureData s = experimentalStructures.get(j);
			//predictedStructures.add(VaryAGamma.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, runEvol, 0.5));
			System.out.println("Good 1 : " + AValue + "\t" + gValue);
			predictedStructures.add(VaryAGamma.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, runEvol, 0.5, AValue, gValue));
			try {
				//Don't want to save full CSV for this time, leave in possibility
				
				//VaryAGamma.saveBenchmarksCSV(new File("output/" + "A_" + AValue + "/results_A_" + AValue+ ".csv"), experimentalStructures, predictedStructures);
				VaryAGamma.saveBenchmarkAvgCSV(avgMetrics,experimentalStructures, predictedStructures, AValue, gValue, j == maxIter - 1);
				j++;
				
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
			BufferedWriter writerAGamma = new BufferedWriter(new FileWriter(outputDirString + "csv/results_A_" + AValue + "_Gamma_" + gValue +evol+ ".csv"));
			for(String str: avgMetrics) {
			  writerAGamma.write(str);
			  String[] line = str.split(",");
			  sensitivity += Double.valueOf(line[2]);
			  ppv += Double.valueOf(line[3]);
			  fscore += Double.valueOf(line[4]);
			  mountainSim += Double.valueOf(line[5]);
			}
			writerAGamma.write("Average,A_" + AValue + "_Gamma_" + gValue +","+ sensitivity/avgMetrics.size()
					+","+ ppv/avgMetrics.size()+","+ fscore/avgMetrics.size()+","+ mountainSim/avgMetrics.size());
			writerAGamma.close();
			
			BufferedWriter writerFinal = new BufferedWriter(new FileWriter(outputDirString + "csv/result"+evol+".csv", true));
			writerFinal.write(AValue +"," + gValue +","+ (sensitivity/avgMetrics.size())+","+ (ppv/avgMetrics.size())+","+ (fscore/avgMetrics.size())+","+ (mountainSim/avgMetrics.size())+"\n");
			writerFinal.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
			
			
	}
	
	public static class AGammaJobInput
	{
		double A;
		double gamma;
		boolean runEvolutionary;
		int uptodataset;
		
		public AGammaJobInput(double a, double g, boolean runEvolutionary,	int uptodataset) {
			this.A = a;
			this.gamma = g;
			this.runEvolutionary = runEvolutionary;
			this.uptodataset = uptodataset;
		}		
	}
	
	class AGammaJobOutput
	{
		
	}	
	
	
    public List<AGammaJobOutput> processInputs(List<AGammaJobInput> inputs)
            throws InterruptedException, ExecutionException {
    	
    	
        int threads = 8;
        ExecutorService service = Executors.newFixedThreadPool(threads);

        List<Future<AGammaJobOutput>> futures = new ArrayList<Future<AGammaJobOutput>>();
        for (final AGammaJobInput input : inputs) {
            Callable<AGammaJobOutput> callable = new Callable<AGammaJobOutput>() {

                public AGammaJobOutput call() throws Exception {
                	AGammaJobOutput output = new AGammaJobOutput();
                	runSpecific(input.A, input.gamma, input.uptodataset, input.runEvolutionary);                   
                    return output;
                }
            };
            futures.add(service.submit(callable));
        }

        service.shutdown();

        List<AGammaJobOutput> outputs = new ArrayList<AGammaJobOutput>();
        for (Future<AGammaJobOutput> future : futures) {
            outputs.add(future.get());
        }
        return outputs;
    }
	
    public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight, double A, double gamma)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");
		
		runEvolutionary = true; 

		// 
		int [] pairedSites = null;
		File structureFile = new File(fastaFile.getAbsolutePath()+".evol.dbn");
		if(!runEvolutionary)
		{
			System.out.println("noevol hit 4");
			structureFile = new File(fastaFile.getAbsolutePath()+".noevol.dbn");
		}
		
		if(structureFile.exists()) // if structure cached, do not bother to recompute
		{
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(structureFile);
			System.out.println("Using cached file " + structureFile);
		}
		else
		{
			System.out.println("Good 2 : " + A + "\t" + gamma);
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
				System.out.println("hit foldOxfold");
				
				new Program().run(argsArray);
				pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".evol.dbn"));		
			}
			else
			{
				String [] argsArray = {fastaFile.getAbsolutePath(), 
						//"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
						"--grammar="+new File("doc/kh.grammar").getAbsolutePath(),
						"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath(),
						"--weight="+weight};			
		
				new Program().run(argsArray);
				System.out.println("noevol hit 5");
				pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
			}	
		}

		StructureData structureData = new StructureData(pairedSites);
		return structureData;		
	}
	
	
	
	public static void saveBenchmarkAvgCSV(ArrayList<String> avgMetrics, List<StructureData> experimentalStructures, List<StructureData> predictedStructures, double AValue, double gamma, boolean isLast) throws IOException
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
			
			avgMetrics.add(AValue +","+ gamma + ","+senstivity+","+ppv+","+fscore+","+mountainSim+predictedStructure.file.getName()+"\n");
			System.out.println(AValue +","+ gamma + ","+senstivity+","+ppv+","+fscore+","+mountainSim+"\n");
		}
		
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
