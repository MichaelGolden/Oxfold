package uk.ac.ox.osscb;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import uk.ac.ox.osscb.analysis.Benchmark;
import uk.ac.ox.osscb.analysis.Benchmarks;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.PPfold;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.visualisation.DataVisualiser;

public class VaryHaasMEA {
	public static void main(String [] args)
	{
		double deltaStart = 1500;
		double deltaEnd = 2500;
		double deltaInc = 100;
		ArrayList<Double> deltaList = new ArrayList<Double>();
		double deltaCurr = deltaStart;
		while(deltaCurr < deltaEnd)
		{
			deltaList.add(deltaCurr);
			deltaCurr += deltaInc;
		}
		
		deltaList.add(Double.POSITIVE_INFINITY);
		File dataDir = new File("datasets/");
		File outputDir = new File("output13/");
		File resultsFile = new File("results_oxfold_evol_test.csv");

		boolean startFile = true;
		
		outputDir.mkdir();
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
		int end = 2;
		File outFile = new File("varyhaas.txt");
		try {
			IO.clearFile(outFile);
		} catch (IOException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		for(double delta : deltaList)
		{
			double sensitivity = 0;
			double ppv = 0;
			double fscore = 0;
			double mountain = 0;
			double t = 0;
			for(int i = 0 ; i < Math.min(end, experimentalStructures.size()) ; i++)
			{
				StructureData s = experimentalStructures.get(i);
				System.out.println(s.file);
				if(datasetno < start)
				{
					datasetno++;
					predictedStructures.add(null);
					continue;
				}
				long startNano = System.nanoTime();
				
			
				StructureData predictedStructure = Benchmarks.foldCofold(outputDir, s.file.getName()+"_cofold", s.sequences, s.sequenceNames, true, delta,Double.POSITIVE_INFINITY, false);
	
				Benchmark benchmark = Benchmark.benchmark(s, predictedStructure);
				System.out.println(benchmark.fscore);
				sensitivity += benchmark.sensitivity;
				ppv += benchmark.ppv;
				fscore += benchmark.fscore;
				mountain += benchmark.mountain;
				t++;
			}
			
			sensitivity /= t;
			ppv /= t;
			fscore /= t;
			mountain /= t;
			 
			try {
				IO.writeLine(outFile, delta+","+sensitivity+","+ppv+","+fscore+","+mountain, true, true);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
	}
}
