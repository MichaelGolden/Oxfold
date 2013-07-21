package uk.ac.ox.osscb.analysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import uk.ac.ox.osscb.CoFoldAnalogue;
import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.KineticFold2;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.Program;

public class VaryKDelta {
	

	
	
	enum FoldingAlgorithm {Oxfold, Cofold, FastOxfold};
	
	public static void main(String [] args)
	{		
		//FoldingAlgorithm algorithm = FoldingAlgorithm.Cofold;
		FoldingAlgorithm algorithm = FoldingAlgorithm.FastOxfold;
		
		boolean runEvol = true;

		int numdatasets = 25;
		if(args.length > 3)
		{
			numdatasets = Integer.parseInt(args[0]);
			double k = Double.parseDouble(args[1]);

			ArrayList<Job> inputs = new ArrayList<Job>();
			for(int i = 2 ; i < args.length ; i++)
			{
				double delta = Double.parseDouble(args[i]);
				//inputs.add(new Job(k, delta, 1, runEvol, numdatasets));
				System.out.println("delta2 "+delta);
				inputs.add(new Job(0.5, k, delta, runEvol, numdatasets)); // for fast oxfold, k=0.5, delta=k, delta2=delta
			}
			VaryKDelta varyKDelta = new VaryKDelta();
			 try {
				List<JobOutput> outputs = varyKDelta.processInputs(inputs, algorithm);
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
			
			/*double [] ks = new double[39];
			double mink = 1 / (double)(ks.length-1);
			for(int i = 0 ; i < ks.length ; i++)
			{
				double alpha = -Math.log10(0.15);
				ks[i] = Math.pow(10, alpha*(2*(mink*i) - 1));
				System.out.print(ks[i]+ " ");
			}
			System.out.println();*/
			
			double [] ks = new double[41];
			double mink = 1 / (double)(ks.length);
			for(int i = 0 ; i < ks.length ; i++)
			{
				double alpha = -Math.log10(0.15);
				ks[i] = (mink)*(i+1)*2000;
				System.out.print(ks[i]+ " ");
			}
			System.out.println();
			
			
			double maxdelta = 2;
			double [] deltas = new double[40];
			double mindelta = 0;
			double incdelta = (maxdelta-mindelta)/((double)deltas.length);
			for(int i = 0 ; i < deltas.length ; i++)
			{
				deltas[i] = mindelta + incdelta*i;
				System.out.print(deltas[i]+ " ");
			}

			System.exit(0);
			
			//double [] ks = {0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
			//double [] deltas = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
			

			ArrayList<Job> inputs = new ArrayList<Job>();
			for(double k : ks)
			{
				for(double delta : deltas)
				{
					inputs.add(new Job(k, delta, 1, runEvol, numdatasets));
				}
			}
			VaryKDelta varyKDelta = new VaryKDelta();
			 try {
				List<JobOutput> outputs = varyKDelta.processInputs(inputs, algorithm);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			/*
			System.err.println("Default: 0 <= K <= 10, 0.1 <= Delta <= 1, with increments of 0.1");
			runAuto(10, 1);*/
		}
		else{
			numdatasets = Integer.parseInt(args[2]);
			Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(args[1])); 
			double delta_o = Constants.IterationCutOff.doubleValue(); 
			double kArg = Double.valueOf(args[0]);
			if(args.length > 2)
			{

				runSpecificOxfold(kArg, delta_o, numdatasets, runEvol);
			}
			else
			{
				runSpecificOxfold(kArg, delta_o, numdatasets, runEvol);
			}
		}
		
	}
	
	private static void runSpecificCofold(double alpha, double tau, int uptodataset, boolean runEvol){
		ArrayList<String> avgMetrics = new ArrayList<String>();
		String outputDirString = "output_cofold/";
		Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(tau)); 
		double delta = Constants.IterationCutOff.doubleValue(); 

		File dataDir = new File("datasets/");
		String evol = runEvol ? "_evol" : "_noevol";
		File outputDir = new File(outputDirString + "K_" + alpha + "_Delta_" + tau +evol+"/");
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
			
			predictedStructures.add(VaryKDelta.foldCofold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, runEvol, alpha, tau, false));
			
			try {
				VaryKDelta.saveBenchmarkAvgCSV(avgMetrics,experimentalStructures, predictedStructures, alpha, delta, j == maxIter - 1);
				j++;
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} //end dataset iterations
		
		try{
			File outputCsvDir = new File(outputDirString + "csv/" );
			outputCsvDir.mkdirs();
			
			double sensitivity = 0; 
			double ppv = 0;
			double fscore = 0; 
			double mountainSim = 0; 
			double n = 0;
			for(int i = 0 ; i < experimentalStructures.size() ; i++)
			{
				String bench = Benchmarks.getBenchmarkString(experimentalStructures.get(i), predictedStructures.get(i));
				String [] split = bench.split(",");
				sensitivity += Double.valueOf(split[2]);
				ppv += Double.valueOf(split[3]);
				fscore += Double.valueOf(split[4]);
				mountainSim += Double.valueOf(split[5]);
				n++;
			}
			
			BufferedWriter writerFinal = new BufferedWriter(new FileWriter(outputDirString + "csv/result"+evol+".csv", true));
			writerFinal.write(alpha +"," + tau +","+ (sensitivity/n)+","+ (ppv/n)+","+ (fscore/n)+","+ (mountainSim/n)+"\n");
			writerFinal.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
			
			
	}
	
	private static void runSpecificOxfold(double KValue, double dValue, int uptodataset, boolean runEvol){
		ArrayList<String> avgMetrics = new ArrayList<String>();
		String outputDirString = "output_oxfold/";
		Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(dValue)); 
		double delta = Constants.IterationCutOff.doubleValue(); 

		File dataDir = new File("datasets/");
		String evol = runEvol ? "_evol" : "_noevol";
		File outputDir = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue +evol+"/");
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
			predictedStructures.add(VaryKDelta.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, runEvol, KValue, 1));
			
			try {
				//Don't want to save full CSV for this time, leave in possibility
				
				
				//VaryKDelta.saveBenchmarksCSV(new File("output/" + "K_" + KValue + "/results_K_" + KValue+ ".csv"), experimentalStructures, predictedStructures);
				VaryKDelta.saveBenchmarkAvgCSV(avgMetrics,experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
				j++;
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} //end dataset iterations
		
		try{
			File outputCsvDir = new File(outputDirString + "csv/" );
			outputCsvDir.mkdirs();
			
			double sensitivity = 0; 
			double ppv = 0;
			double fscore = 0; 
			double mountainSim = 0; 
			double n = 0;
			for(int i = 0 ; i < experimentalStructures.size() ; i++)
			{
				String bench = Benchmarks.getBenchmarkString(experimentalStructures.get(i), predictedStructures.get(i));
				String [] split = bench.split(",");
				sensitivity += Double.valueOf(split[2]);
				ppv += Double.valueOf(split[3]);
				fscore += Double.valueOf(split[4]);
				mountainSim += Double.valueOf(split[5]);
				n++;
			}
			
			BufferedWriter writerFinal = new BufferedWriter(new FileWriter(outputDirString + "csv/result"+evol+".csv", true));
			writerFinal.write(KValue +"," + dValue +","+ (sensitivity/n)+","+ (ppv/n)+","+ (fscore/n)+","+ (mountainSim/n)+"\n");
			writerFinal.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
			
			
	}
	
	static Random random= new Random(983118011408130222L);
	/*public static void selectN(StructureData s, int n)
	{		
		int select = Math.min(s.sequences.size(), n);
		while(select < s.sequences.size())
		{
			int i = random.nextInt(s.sequences.size());
			s.sequences.remove(i);
			s.sequenceNames.remove(i);
		}
		
		boolean [] delete = new boolean[s.sequences.get(0).length()];
		Arrays.fill(delete, true);
		for(int i = 0 ; i < s.sequences.size() ; i++)
		{
			String seq = s.sequences.get(i);
			for(int j = 0 ; j < seq.length() ; j++)
			{
				if(seq.charAt(j) != '-')
				{
					delete[j] = false;
				}
			}
		}
		
		for(int i = 0 ; i < delete.length ; i++)
		{
			if(delete[i])
			{
				int pairedTo = s.pairedSites[i]-1;
				if(pairedTo >= 0)
				{
					s.pairedSites[i] = 0;
					s.pairedSites[pairedTo]=0;
				}				
			}
		}
	}*/
	
	public static void selectN(ArrayList<String> sequences, ArrayList<String> sequenceNames, int n)
	{
		
		int select = Math.min(sequences.size(), n);
		while(select < sequences.size())
		{
			int i = random.nextInt(sequences.size());
			sequences.remove(i);
			sequenceNames.remove(i);
		}
	}
	

    int threads = 8;
	private static void runSpecificFastOxfold(double KValue, double dValue, double delta2, int uptodataset, boolean runEvol){
		ArrayList<String> avgMetrics = new ArrayList<String>();
		String outputDirString = "output_rrna/";
		
		Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(dValue)); 
		
		double delta = Constants.IterationCutOff.doubleValue(); 

		File dataDir = new File("datasets/");
		String evol = runEvol ? "_evol" : "_noevol";
		File outputDir = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue+"_delta2_"+delta2 +evol+"/");
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
		sublist.addAll((List<StructureData>)experimentalStructures.subList(0, Math.min(experimentalStructures.size(), maxIter)));
		experimentalStructures = sublist;
		
		ArrayList<StructureData> predictedStructures = new ArrayList<StructureData>();
		ArrayList<StructureData> noDelta2 = new ArrayList<StructureData>();
		 
		ArrayList<StructureData> usedExperimental = new ArrayList<StructureData>();
		int skip = 0;
		int j = 0;
		for(int i = 0 ; i < experimentalStructures.size() ; i++) 
		{
			if(i < skip)
			{
				continue;
			}
			System.out.println(experimentalStructures.get(i).file);
			StructureData s= experimentalStructures.get(i);
			selectN(s.sequences, s.sequenceNames, (int)delta2);
			System.out.println("VAL"+s.sequences.size()+"\t"+delta2);
			usedExperimental.add(s);
			
			predictedStructures.add(foldOxfold(outputDir, s.file.getName()+"_", s.sequences, s.sequenceNames, runEvol, KValue, delta2));
			//noDelta2.add(foldOxfold(outputDir, s.file.getName()+"_", s.sequences, s.sequenceNames, runEvol, KValue, delta2));
			//File dir2 = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue+"_delta2_"+(1.0) +evol+"/");
			//dir2.mkdirs();
			//noDelta2.add(foldOxfold(dir2, s.file.getName(), s.sequences, s.sequenceNames, runEvol, KValue, 1));
			
			try {				
				VaryKDelta.saveBenchmarkAvgCSV(avgMetrics,experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
				j++;
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} //end dataset iterations
		
		try{
			File outputCsvDir = new File(outputDirString + "csv/" );
			outputCsvDir.mkdirs();
			
			double sensitivity = 0; 
			double ppv = 0;
			double fscore = 0; 
			double mountainSim = 0; 
			double n = 0;
			for(int i = 0 ; i < usedExperimental.size() ; i++)
			{
				String bench = Benchmarks.getBenchmarkString(usedExperimental.get(i), predictedStructures.get(i));
				String [] split = bench.split(",");
				sensitivity += Double.valueOf(split[2]);
				ppv += Double.valueOf(split[3]);
				fscore += Double.valueOf(split[4]);
				mountainSim += Double.valueOf(split[5]);
				n++;
			}
			

			double iterWithDelta2 = 0;
			double iterWithoutDelta2 = 0;
			for(int i = 0 ; i < predictedStructures.size() ; i++)
			{
				//System.out.println("ACVB "+predictedStructures.get(i).title);
				String with = predictedStructures.get(i).title.split(";")[1].split("=")[1].split(",")[0];
				//System.out.println(">>>"+i+"\t"+with+"\t"+without);
				iterWithDelta2 += Double.parseDouble(with); 
				

				if(i < noDelta2.size())
				{
				String without = noDelta2.get(i).title.split(";")[1].split("=")[1].split(",")[0];
				iterWithoutDelta2 += Double.parseDouble(without);
				}
			}
			
			
			
			BufferedWriter writerFinal = new BufferedWriter(new FileWriter(outputDirString + "csv/result"+evol+".csv", true));
			writerFinal.write(KValue +"," + dValue+","+delta2 +","+ (sensitivity/n)+","+ (ppv/n)+","+ (fscore/n)+","+ (mountainSim/n)+"\t"+(iterWithDelta2/n)+"\t"+(iterWithoutDelta2/n)+"\n");
			writerFinal.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			
			
			
	}
	
	public static class Job
	{
		double v1;
		double v2;
		double v3;
		boolean runEvolutionary;
		int uptodataset;
		
		public Job(double k, double delta, double delta2, boolean runEvolutionary,	int uptodataset) {
			this.v1 = k;
			this.v2 = delta;
			this.v3 = delta2;
			this.runEvolutionary = runEvolutionary;
			this.uptodataset = uptodataset;
		}		
	}
	
	class JobOutput
	{
		
	}	
	
    public List<JobOutput> processInputs(List<Job> inputs, final FoldingAlgorithm algorithm)
            throws InterruptedException, ExecutionException {

        ExecutorService service = Executors.newFixedThreadPool(threads);

        List<Future<JobOutput>> futures = new ArrayList<Future<JobOutput>>();
        for (final Job input : inputs) {
            Callable<JobOutput> callable = new Callable<JobOutput>() {

                public JobOutput call() throws Exception {
                	JobOutput output = new JobOutput();
                	switch(algorithm)
                	{
                		case FastOxfold:
                			runSpecificFastOxfold(input.v1, input.v2, input.v3, input.uptodataset, input.runEvolutionary);
                		break;
	                	case Oxfold:
	                		runSpecificOxfold(input.v1, input.v2, input.uptodataset, input.runEvolutionary);
	                		break;
	                	case Cofold:
	                		runSpecificCofold(input.v2, input.v1, input.uptodataset, input.runEvolutionary);
	                		break;
                	}
                    return output;
                }
            };
            futures.add(service.submit(callable));
        }

        service.shutdown();

        List<JobOutput> outputs = new ArrayList<JobOutput>();
        for (Future<JobOutput> future : futures) {
            outputs.add(future.get());
        }
        return outputs;
    }
	
	private static  void runAuto(double k, double delta){
		ArrayList<String> avgMetrics = new ArrayList<String>();
		String outputDirString = "output/";
		//variables for K or delta and their ranges. Then move below out to new method. 
				for(double KValue = 0; KValue <= k; KValue += 0.1){ //change this later
					for(double dValue = 0; dValue <= delta; dValue += 0.1){
						Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(dValue)); 
						delta = Constants.IterationCutOff.doubleValue(); 
						if(KValue == 0){
							KValue = 0.1; 
						}
						File dataDir = new File("datasets/");
						File outputDir = new File(outputDirString + "K_" + KValue + "_Delta_" + dValue +"/");
						outputDir.mkdir();
						
						ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
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
						

						ArrayList<StructureData> predictedStructures = new ArrayList<StructureData>();
						
						int maxIter = 1; 
						for(int j  = 0; j < maxIter; j++) //go through select # of datasets
						//for(StructureData s : experimentalStructures) //go through all
						{
							StructureData s = experimentalStructures.get(0);
							//StructureData s = experimentalStructures.get(j);
							predictedStructures.add(VaryKDelta.foldOxfold(outputDir, s.file.getName(), s.sequences, s.sequenceNames, false, KValue, 1));
							
							try {
								//Don't want to save full CSV for this time, leave in possibility
								
								
								//VaryKDelta.saveBenchmarksCSV(new File("output/" + "K_" + KValue + "/results_K_" + KValue+ ".csv"), experimentalStructures, predictedStructures);
								VaryKDelta.saveBenchmarkAvgCSV(avgMetrics, experimentalStructures, predictedStructures, KValue, delta, j == maxIter - 1);
								
								
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
	
	public static StructureData foldCofold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double alpha, double tau, boolean reverseGrammar)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");
		File basePairProbFile = null;

		int [] pairedSites = null;
		File structureFile = new File(fastaFile.getAbsolutePath()+".evol.dbn");
		if(!runEvolutionary)
		{
			structureFile = new File(fastaFile.getAbsolutePath()+".noevol.dbn");
		}
		
		if(structureFile.exists()) // if structure cached, do not bother to recompute
		{
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(structureFile);
			System.out.println("Using cached file " + structureFile);
		}
		else
		{
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
	
			
		}
		StructureData structureData = new StructureData(pairedSites);
		if(basePairProbFile != null && basePairProbFile.exists())
		{
			structureData.basePairProbFile =  basePairProbFile;
		}
		return structureData;		
	}
	
	
	public static StructureData foldOxfold(File dir, String name, List<String> sequences, List<String> sequenceNames, boolean runEvolutionary, double weight, double delta2)
	{
		File fastaFile = new File(dir.getAbsolutePath()+File.separatorChar+name+".fas");
		IO.saveToFASTAfile(sequences, sequenceNames, fastaFile);
		File outNewick = new File(dir.getAbsolutePath()+File.separatorChar+name+".nwk");

		// 
		int [] pairedSites = null;
		File structureFile = new File(fastaFile.getAbsolutePath()+".evol.dbn");
		if(!runEvolutionary)
		{
			structureFile = new File(fastaFile.getAbsolutePath()+".noevol.dbn");
		}
		
		double elapsedNano = -1;
		
		String title = "";
		if(structureFile.exists()) // if structure cached, do not bother to recompute
		{
			title = RNAFoldingTools.getLine(structureFile, 0);
			pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(structureFile);
			System.out.println("Using cached file " + structureFile);
		}
		else
		{
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
		
				//new Program().run(argsArray);
				long startNano = System.nanoTime();
				new KineticFold2().foldEvolutionary(fastaFile.getAbsolutePath(),new File("doc/ppfold.grammar").getAbsolutePath(),new File("doc/ppfold.parameters").getAbsolutePath(),outNewick.getAbsolutePath(),weight, delta2);
				long endNano = System.nanoTime();
				elapsedNano = (endNano - startNano)/1000000000.0;
				
				title = RNAFoldingTools.getLine(structureFile, 0);
				pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(structureFile);		
			}
			else
			{
				String [] argsArray = {fastaFile.getAbsolutePath(), 
						//"--grammar="+new File("doc/kh_reverse.grammar").getAbsolutePath(),
						"--grammar="+new File("doc/kh.grammar").getAbsolutePath(),
						"--grammar-params="+new File("doc/kh.parameters").getAbsolutePath(),
						"--weight="+weight};			
		
				new Program().run(argsArray);
				
				title = RNAFoldingTools.getLine(structureFile, 0);
				pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketFile(new File(fastaFile.getAbsolutePath()+".noevol.dbn"));		
			}	
		}

		StructureData structureData = new StructureData(pairedSites);
		structureData.time = elapsedNano;
		structureData.title=  title;
		return structureData;		
	}
	
	
	public static void saveBenchmarkAvgCSV(ArrayList<String> avgMetrics, List<StructureData> experimentalStructures, List<StructureData> predictedStructures, double KValue, double delta, boolean isLast) throws IOException
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
