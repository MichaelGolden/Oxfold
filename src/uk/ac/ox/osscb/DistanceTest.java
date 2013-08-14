package uk.ac.ox.osscb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

import uk.ac.ox.osscb.analysis.AlignmentUtils;
import uk.ac.ox.osscb.analysis.Benchmark;
import uk.ac.ox.osscb.analysis.Benchmarks;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.PPfold;
import uk.ac.ox.osscb.analysis.RankingAnalyses;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.analysis.VaryKDelta;

public class DistanceTest {
	public static void main(String[] args) {

		File dataDir = new File("datasets/");
		File outputDir = new File("output_long/");
		//File benchmarksFile = new File("benchmarks.csv");
		//File summaryFile = new File("summary.csv");
		/*try {
			IO.clearFile(benchmarksFile);
			IO.clearFile(summaryFile);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}*/


		outputDir.mkdir();
		System.out.println(dataDir.list().length);
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for (File experimentalFile : dataDir.listFiles()) {
			try {
				StructureData data = StructureData
						.readExperimentalStructureData(experimentalFile);
				if (data.pairedSites.length > 100
						&& experimentalFile.getName().contains("TestRNAData")) {
					continue;
				}
				experimentalStructures.add(data);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		Collections.sort(experimentalStructures);

		int m = 0;
		for (StructureData e : experimentalStructures) {
			System.out.println(m + "\t" + e.file + "\t"
					+ e.sequences.get(0).length());
			m++;
		}

		ArrayList<StructureData> predictedStructures = new ArrayList<StructureData>();
		int datasetno = 0;
		// int start = 0;
		// int end = 1;
		int start = 0;
		int end = 1000;

		double [] taus = {25,600, 700,500,800,400,900,300,1000,200,1200,100,1400,50, 1700,2000,2500,150,1300};
		
		for(double tau : taus)
		{
			for (int i = 0; i < Math.min(end, experimentalStructures.size()); i++) {
				StructureData s = experimentalStructures.get(i);
				if (datasetno < start) {
					datasetno++;
					predictedStructures.add(null);
					continue;
				}
	
	
				File dir = new File(outputDir + "/" + s.file.getName());
				dir.mkdirs();
				double alpha = 0;
				StructureData ppfold = VaryKDelta.foldCofold(dir,  s.file.getName()+"ppfold_analogue", s.sequences, s.sequenceNames, true, alpha, Double.POSITIVE_INFINITY, false);

				StructureData oxfold = VaryKDelta.foldCofold(dir,  s.file.getName()+"_tau_"+tau, s.sequences, s.sequenceNames, true, alpha, tau, false);
				try {
					File outFile = new File(outputDir+"/distance_"+tau+".csv"); 
					System.out.println(outFile);
					IO.writeLine(outFile, Benchmarks.getBenchmarkString(s, oxfold)+","+Benchmarks.getBenchmarkString(s, ppfold), true, true);
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
		}
	}

}
