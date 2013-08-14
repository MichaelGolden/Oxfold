package uk.ac.ox.osscb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import uk.ac.ox.osscb.analysis.AlignmentUtils;
import uk.ac.ox.osscb.analysis.Benchmark;
import uk.ac.ox.osscb.analysis.Benchmarks;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.PPfold;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.RankingAnalyses;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.analysis.VaryKDelta;
import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.visualisation.DataVisualiser;
import uk.ac.ox.osscb.visualisation.SVG;

public class SimonsIdea {
	public static void main(String[] args) {
		
		File dataDir = new File("large_more/");
		File outputDir = new File("output_test2/");
		File resultsFile = new File("results_oxfold_evol_test.csv");
		File benchmarksFile = new File("benchmarks.csv");
		File summaryFile = new File("summary.csv");
		try {
			IO.clearFile(benchmarksFile);
			IO.clearFile(summaryFile);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		boolean startFile = true;
		
		outputDir.mkdir();
		System.out.println(dataDir.list().length);
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for(File experimentalFile : dataDir.listFiles())
		{
			try {
				StructureData data = StructureData.readExperimentalStructureData(experimentalFile);
				//System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(data.pairedSites));
				//System.out.println(getHelices(data.pairedSites));
				//System.exit(0);
				if(data.pairedSites.length > 100 && experimentalFile.getName().contains("TestRNAData"))
				{
					continue;
				}
				
				/*while(data.sequences.size() > 1)
				{
					data.sequences.remove(data.sequences.size()-1);
					data.sequenceNames.remove(data.sequenceNames.size()-1);
				}
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
		int end = 40;
		DecimalFormat df = new DecimalFormat("0.0000");
		/*File outputFile = new File("test.txt");
		try {
			IO.clearFile(outputFile);
		} catch (IOException e3) {
			// TODO Auto-generated catch block
			e3.printStackTrace();
		}
		
	
		double maxdelta2 = 1;
		double mindelta2 = 0;
		double incdelta2 = 0.025;
		double current2 = mindelta2;
		ArrayList<Double> alphas = new ArrayList<Double>();
		while(current2 < maxdelta2)
		{
			alphas.add(current2);
			current2 += incdelta2;
		}
		*/
		
		//for(double alpha : alphas)
		//{
		double total = 0;
		double n = 0;
		
		double maxdelta = 1;
		double mindelta = -0.95;
		double incdelta = 0.05;
		double current = mindelta;
		ArrayList<Double> deltas = new ArrayList<Double>();
		while(current < maxdelta)
		{
			deltas.add(current);
			current += incdelta;
		}
		
		ArrayList<BenchmarkStore> stores = new ArrayList<BenchmarkStore>(); 
		for(int i = 0 ; i < Math.min(end, experimentalStructures.size()) ; i++)
		{
			StructureData s = experimentalStructures.get(i);
			if(datasetno < start)
			{
				datasetno++;
				predictedStructures.add(null);
				continue;
			}
			
			//StructureData predictedStructure = Benchmarks.foldOxfold(outputDir, s.file.getName()+"_oxfold", s.sequences, s.sequenceNames, true,0.5, false);
			//StructureData predictedStructure2 = Benchmarks.foldRNAalifold(outputDir, s.file.getName()+"_rnaalifold", s.sequences, s.sequenceNames);
			//StructureData predictedStructure3 = PPfold.fold(outputDir, s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true);

			
			
			ArrayList<Benchmark> benchmarks = new ArrayList<Benchmark>();

			File dir = new File(outputDir+"/"+s.file.getName());
			dir.mkdirs();
			StructureData ppfold = PPfold.fold(dir,  s.file.getName()+"_ppfold", s.sequences, s.sequenceNames, true, true);
			Benchmark benchmarkPPfold = Benchmark.benchmark(s, ppfold);
			boolean first = true;
			ArrayList<Double> meaDeltas =new ArrayList<Double>();
			ArrayList<Double> matrixDeltas = new ArrayList<Double>();
			double avgDiversity = AlignmentUtils.calculateNucleotideDiversityAverage(s.sequences);
			double reliability = 0;

			double [][] deltaMatrix = null;
			BenchmarkStore currentBenchmark = new BenchmarkStore();
			for(double delta : deltas)
			{					
				double cdelta = Double.parseDouble(df.format(delta));
				Constants.IterationCutOff = PointRes.valueOf(Double.valueOf(cdelta));
				Constants.IterationCutOffDouble = cdelta;
				
			
				System.out.println(s.file.getName());
				StructureData oxfold = null;

				try
				{
					oxfold = VaryKDelta.foldOxfold(dir,  s.file.getName()+"_delta_"+df.format(cdelta), s.sequences, s.sequenceNames, true, Double.POSITIVE_INFINITY, 0);
				}
				catch(Exception ex)
				{
					ex.printStackTrace();
				}
				

				if(first)
				{
					try {
						BufferedReader buffer = new BufferedReader(new FileReader(outputDir+"/"+ s.file.getName()+"/"+ s.file.getName()+"_delta_"+df.format(cdelta)+".fas.deltas"));
						String structure = buffer.readLine();
						double [][] bp = StructureData.getBasePairProb(ppfold.basePairProbFile);
						double [] single = StructureData.getSingleProbs(bp);
						deltaMatrix = PosteriorProbabilitiesCalculator.getDiffs(bp);
						for(int x = 0 ; x < deltaMatrix.length ; x++)
						{
							for(int y = 0 ; y < deltaMatrix.length ; y++)
							{
								System.out.print(deltaMatrix[x][y]+"\t");
							}
							System.out.println();
						}

						for(int k = 0 ; k < ppfold.pairedSites.length ; k++)
						{
							if(ppfold.pairedSites[k] == 0)
							{
								reliability += single[k];
							}
							else
							{
								reliability += bp[k][ppfold.pairedSites[k]-1];
							}
						}
						reliability /= (double)(ppfold.pairedSites.length);
						
						String [] split = buffer.readLine().split(",");
						for(String text : split)
						{
							meaDeltas.add(Double.parseDouble(text));
						}

						/*
						double[][] matrix = new double[oxfold.pairedSites.length][oxfold.pairedSites.length];

						
						String textline = null;
						for(int x = 0 ; x < matrix.length  && (textline = buffer.readLine()) != null ; x++)
						{
							String [] split2 = textline.split(",");
							for(int y = 0 ; y < split2.length ; y++)
							{
								double value = Double.parseDouble(split2[y]);
								if(value != -1)
								{
									matrix[x][y] = value;			
									matrixDeltas.add(value);
								}
							}
						}*/
						for(int x = 0 ; x < bp.length  ; x++)
						{
							for(int y = 0 ; y < bp.length ; y++)
							{
								if(bp[x][y] > 0)
								{
									matrixDeltas.add(bp[x][y]);
								}
							}
						}
						
						buffer.close();
						
					} catch (IOException e2) {
						// TODO Auto-generated catch block
						e2.printStackTrace();
					}
					first=false;
				}
				
				if(oxfold != null)
				{
					Benchmark benchmarkOxfold = Benchmark.benchmark(s, oxfold);
					benchmarkOxfold.param1 = delta;
					benchmarks.add(benchmarkOxfold);
					currentBenchmark.add(delta, benchmarkOxfold.fscore);
					currentBenchmark.reliabilityScore = reliability;
					currentBenchmark.normalisedEntropy = ppfold.normalisedEntropy;
					stores.add(currentBenchmark);
					try {
						IO.writeLine(benchmarksFile, cdelta+","+s.sequences.size()+","+benchmarkOxfold.toString()+","+benchmarkPPfold.toString(), true, true);
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
				}
			}
			
			double maxFscore = 0;
			int startpos = 0;
			for(int a = 0 ; a < benchmarks.size() ; a++)
			{
				if(benchmarks.get(a).fscore > maxFscore)
				{
					maxFscore = benchmarks.get(a).fscore;
					startpos = a;
				}
			}
			int endpos = startpos;
			/*double scale = 0.95;
			while(startpos-1 >= 0 && benchmarks.get(startpos-1).fscore >= maxFscore*scale)
			{
				startpos--;
			}*/
			for(int a = startpos ; a < benchmarks.size() ; a++)
			{
				if(benchmarks.get(a).fscore == maxFscore)
				//if(benchmarks.get(a).fscore >= maxFscore*scale)
				{
					endpos = a;
				}
			}
			int middle = (startpos+endpos)/2;
			try {
				String structureDist = "";
				String matrixDist = "";
				try
				{
					structureDist = RankingAnalyses.getMin(meaDeltas)+","+RankingAnalyses.getPercentile(meaDeltas, 0.125)+","+RankingAnalyses.getMedian(meaDeltas)+","+RankingAnalyses.getPercentile(meaDeltas, 90);
					matrixDist = RankingAnalyses.getMin(matrixDeltas)+","+RankingAnalyses.getPercentile(matrixDeltas, 0.2)+","+RankingAnalyses.getMedian(matrixDeltas)+","+RankingAnalyses.getPercentile(matrixDeltas, 0.8)+","+RankingAnalyses.getPercentile(matrixDeltas, 0.9)+","+RankingAnalyses.getMax(matrixDeltas)+","+reliability;
				}
				catch(Exception ex)
				{
					structureDist = ","+","+",";
					matrixDist = ","+","+","+","+","+",";
					ex.printStackTrace();
				}
				
				double minDelta = Double.MAX_VALUE;
				ArrayList<Helix> helices = getHelices(ppfold.pairedSites, deltaMatrix);
				for(Helix helix : helices)
				{
					minDelta = Math.min(minDelta, helix.maxdelta);
				}

				currentBenchmark.minDelta = minDelta;
				
				IO.writeLine(summaryFile,deltas.get(startpos)+","+deltas.get(middle)+","+ deltas.get(endpos)+","+s.sequences.size()+","+ppfold.entropy+","+ppfold.normalisedEntropy+","+reliability+","+minDelta+","+structureDist+","+matrixDist+","+benchmarks.get(middle).toString()+","+benchmarkPPfold.toString(), true, true);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			/*
			double pick = RankingAnalyses.getPercentile(meaDeltas, alpha);
			double smallest = Double.MAX_VALUE;
			int chosen = 0;
			for(int k = 0 ; k < deltas.size() ; k++)
			{
				if(Math.abs(deltas.get(k) - pick) < smallest)
				{
					smallest = Math.abs(deltas.get(k) - pick);
					chosen = k;
				}
			}
			
			try {
				total += benchmarks.get(chosen).fscore;
				n++;
				IO.writeLine(outputFile, chosen+","+benchmarks.get(chosen).fscore+",", false, true);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			*/
			
			/*
			try
			{
				BufferedWriter writer = new BufferedWriter(new FileWriter(benchmarksFile, true));
				for(Benchmark benchmark : benchmarks)
				{
					writer.write(benchmark.toString()+"\n");
				}
				writer.close();
			}
			catch(IOException ex)
			{
				ex.printStackTrace();
			}*/
		}
		/*try {
			IO.writeLine(outputFile, "Alpha,"+alpha+ ",Average,"+(total/n), true, true);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		}*/		
		
		
		
		for(double delta : deltas)
		{
			double sum = 0;
			double count = 0;
			
			for(BenchmarkStore s : stores)
			{
			
	
				sum += s.pickBestFscore(delta);
				count++;
			}
		
			double average = sum/count;
			System.out.println("D\t"+delta+"\t"+average);
		}
		
		Random random = new Random(139194111433014444L);
		double max = 0;
		while(true)
		{
			double a =  random.nextDouble()*2-1;
			double b =  random.nextDouble()*2-1;
			double c =  random.nextDouble()*2-1;
			
			double sum = 0;
			double count = 0;
			for(BenchmarkStore s : stores)
			{
			
				double calc = a*s.minDelta+b*s.reliabilityScore+c;
				double delta = Math.min(1,Math.max(-1, calc));
	
				sum += s.pickBestFscore(delta);
				count++;
			}
			double average = sum/count;
			if(average > max)
			{
				max = average;
				System.out.println(a+"\t"+b+"\t"+c+"\t"+average);
			}
		}
		
	}
	
	public static class BenchmarkStore
	{
		double minDelta;
		double normalisedEntropy;
		double reliabilityScore;
		int n;
		
		ArrayList<Double> deltas = new ArrayList<Double>();
		ArrayList<Double> fscores = new ArrayList<Double>();

		public void add(double delta, double fscore)
		{
			this.deltas.add(delta);
			this.fscores.add(fscore);
		}
		
		public double pickBestFscore(double delta)
		{
			int pick = pickNearestDelta(delta);
			return fscores.get(pick);
		}
		
		public int pickNearestDelta(double delta)
		{
			double smallest = Double.MAX_VALUE;
			int chosen = 0;
			for(int k = 0 ; k < deltas.size() ; k++)
			{
				if(Math.abs(deltas.get(k) - delta) < smallest)
				{
					smallest = Math.abs(deltas.get(k) - delta);
					chosen = k;
				}
			}
			
			return chosen;
		}
	}
	
	public static ArrayList<Helix> getHelices(int [] pairedSites, double [][] deltas)
	{
		ArrayList<Helix> helices =  getHelices(pairedSites);
		
		for(Helix helix : helices)
		{
			double maxdelta = 0;
			for(int i = 0 ; i < helix.getHelixLength() ; i++)
			{
				maxdelta = Math.max(maxdelta, deltas[helix.getLeftIdx()+i][helix.getRightIdx()-i]);
			}
			helix.maxdelta = maxdelta;
			System.out.println("HELIX"+helix+"\t"+maxdelta);
		}
		
		return helices;
	}
	
	
	public static ArrayList<Helix> getHelices(int [] pairedSites)
	{
		ArrayList<Helix> helices = new ArrayList<Helix>();
		for(int i = 0 ; i < pairedSites.length ; i++)
		{
			int x = i;
			int y = pairedSites[i];
			int length = 0;
			
			if(y > x)
			{
				int j = i;
				for(; j < pairedSites.length ; j++)
				{
					if(pairedSites[j] == 0)
					{
						break;
					}
					else
					{
						length++;
					}
				}
				i = j;
				
				helices.add(new Helix(x,y-1,length));
			}
		}
		
		return helices;
		
	}

}
