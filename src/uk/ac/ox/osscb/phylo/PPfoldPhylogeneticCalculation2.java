package uk.ac.ox.osscb.phylo;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import uk.ac.ox.osscb.analysis.AmbiguityCodes;
import uk.ac.ox.osscb.analysis.IO;

/**
 * After reading and parsing input, algorithm execution starts here. (Note:
 * Columns with many gaps are already removed)
 * 
 * @author Z.Sukosd
 */
public class PPfoldPhylogeneticCalculation2 {
	
	public static void main(String [] args)
	{
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File("TestRNAData14.dat_cofold.fas"), sequences, sequenceNames);
		try {
			double [][] probs = getPhyloProbs(sequences, sequenceNames, new File("TestRNAData14.dat_cofold.nwk"));
			/*for(int i = 0 ; i < probs.length ; i++)
			{
				for(int j = 0 ; j < probs.length ; j++)
				{
					System.out.print(probs[i][j]+"\t");
				}
				System.out.println();
			}*/
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static double [][] getPhyloProbs(ArrayList<String> sequences, ArrayList<String> sequenceNames, File newickFile) throws Exception
	{
		ArrayList<char[]>  columns = new ArrayList<char[]>();
		int columnCount = sequences.get(0).length();
		for(int j = 0 ; j < columnCount ; j++)
		{
			char [] column = new char[sequences.size()];
			columns.add(column);
		}
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			String seq = sequences.get(i);
			for(int j = 0 ; j < seq.length() ; j++)
			{
				columns.get(j)[i] += seq.charAt(j);  
			}
			
		}
		
		
		/*for(char[] column : columns)
		{
			for(int i = 0 ; i < column.length ; i++)
			{
				System.out.print(column[i]);
			}
			System.out.println();
		}
		*/
		String parametersFile = "doc/matrices.in";
		Parameters param = null;
		try
		{
			param = Parameters.readParam(parametersFile);
		}
		catch(Exception ex)
		{
			param = Parameters.cached.get(parametersFile);
			if(param == null)
			{
				throw ex;
			}
		}
		AsynchronousJobExecutor executor = new AsynchronousJobExecutorBlocking();
		Tree  tree = NewickReader.readNewick(newickFile.getAbsolutePath());
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		
		//long starttime = System.currentTimeMillis();
		double[][] probmatrix = PPfoldPhylogeneticCalculation2.createPhyloProb(NullProgress.INSTANCE
				.getChildProgress(phylopart), 1, tree, columns,
				sequenceNames, columns.size(), param, executor, false, 0);
		/*System.out.println("probmatrix!!!");
		for(int i = 0; i < probmatrix.length; i++){
			for (int j = 0; j < probmatrix[i].length; j++){
				System.out.print(probmatrix[i][j] + "\t");
			}
			System.out.println(); 
		}*/
		return probmatrix;
	}
	
	public static double [][] getPhyloProbs(ArrayList<String> sequences, ArrayList<String> sequenceNames, File newickFile, double p) throws Exception
	{
		ArrayList<char[]>  columns = new ArrayList<char[]>();
		int columnCount = sequences.get(0).length();
		for(int j = 0 ; j < columnCount ; j++)
		{
			char [] column = new char[sequences.size()];
			columns.add(column);
		}
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			String seq = sequences.get(i);
			for(int j = 0 ; j < seq.length() ; j++)
			{
				columns.get(j)[i] += seq.charAt(j);  
			}
			
		}
		String parametersFile = "doc/matrices.in";
		Parameters param = null;
		try
		{
			param = Parameters.readParam(parametersFile);
		}
		catch(Exception ex)
		{
			param = Parameters.cached.get(parametersFile);
			if(param == null)
			{
				throw ex;
			}
		}
		AsynchronousJobExecutor executor = new AsynchronousJobExecutorBlocking();
		Tree  tree = NewickReader.readNewick(newickFile.getAbsolutePath());
		
		double scfg_to_phylo_ratio = 0.095 * columns.size()
				/ columns.get(0).length;
		double phylopart = 0.95 / (scfg_to_phylo_ratio + 1);
		double scfgpart = 0.95 - phylopart; // max is 95% because processing
		double[][] probmatrix = PPfoldPhylogeneticCalculation2.createPhyloProb(NullProgress.INSTANCE
				.getChildProgress(phylopart), 1, tree, columns,
				sequenceNames, columns.size(), param, executor, false, 0);
		AmbiguityCodes decoder = new AmbiguityCodes();
		
		int tmpCount = 0; 
		
		for(int a = 0; a < columns.size(); a++)
		{
			for(int b = 0; b< columns.size(); b++){
				HashSet<Character> diffCol1 = new HashSet<Character>(); 
				HashSet<Character> diffCol2 = new HashSet<Character>(); 
				/*
				HashMap<Character, Double> diffColMap1 = new HashMap<Character, Double>(); 
				HashMap<Character, Double> diffColMap2 = new HashMap<Character, Double>(); 
				int countAdded1 = 0; 
				int countAdded2 = 0; 
				
				for(int i = 0 ; i < columns.get(a).length ; i++){
					char [] toAdd = decoder.getAmbigChars(columns.get(a)[i]);
					countAdded1++; 
					for(int ii = 0; ii < toAdd.length; ii++){
						if(diffColMap1.containsKey(toAdd[ii])){
							diffColMap1.put(toAdd[ii], diffColMap1.get(toAdd[ii]) + (1.0 / toAdd.length));
							System.out.println((1.0 / toAdd.length));
						}
						else{
							diffColMap1.put(toAdd[ii], (1.0 / toAdd.length));
						}
					}
				}
				for(int i = 0 ; i < columns.get(b).length ; i++){
					char [] toAdd = decoder.getAmbigChars(columns.get(a)[i]);
					countAdded2++; 
					for(int ii = 0; ii < toAdd.length; ii++){
						if(diffColMap2.containsKey(toAdd[ii])){
							diffColMap2.put(toAdd[ii], diffColMap2.get(toAdd[ii]) + (1.0 / toAdd.length));
							System.out.println((1.0 / toAdd.length));
						}
						else{
							diffColMap2.put(toAdd[ii], (1.0 / toAdd.length));
						}
					}
				}
				*/
				
				//Decode ambiguous characters down to the normal 4 RNA bases then add into set
				for(int i = 0 ; i < columns.get(a).length ; i++)
				{
					char [] toAdd = decoder.getAmbigChars(columns.get(a)[i]);
					for(int ii = 0; ii < toAdd.length; ii++){
						diffCol1.add(toAdd[ii]);
					}
				}
				for(int i = 0 ; i < columns.get(b).length ; i++)
				{
					char [] toAdd = decoder.getAmbigChars(columns.get(b)[i]);
					for(int ii = 0; ii < toAdd.length; ii++){
						diffCol2.add(toAdd[ii]);
					}
				}
				
				
				/*
				//Original
				for(int i = 0 ; i < columns.get(a).length ; i++)
				{
					diffCol1.add(columns.get(a)[i]);
				}
				for(int i = 0 ; i < columns.get(b).length ; i++)
				{
					diffCol2.add(columns.get(b)[i]);
				}
				*/
				
				if((diffCol1.size() != diffCol2.size()) && (diffCol1.size() == 1 || diffCol2.size() == 1 )){
					probmatrix[a][b] = Math.max(0, probmatrix[a][b] * p );//- p);//penalty for certain paired with uncertain
					tmpCount++;
				}
			}
		}
		/*System.out.println("penalty hits = " + tmpCount);
		System.out.println("probmatrix!!!");
		for(int i = 0; i < probmatrix.length; i++){
			for (int j = 0; j < probmatrix[i].length; j++){
				System.out.print(probmatrix[i][j] + "\t");
			}
			System.out.println(); 
		}
		

		System.out.println("GOT P THROUGH! " +p);
		*/
		return probmatrix;
	}

	private static double[][] createPhyloProb(Progress act, int userjobsnr,
			Tree tree, List<char[]> columns_char, List<String> names,
			final int length, Parameters param,
			AsynchronousJobExecutor executor, final boolean verbose, int execnr)
			throws InterruptedException {
		final long starttime = System.nanoTime();
		if(verbose){
		System.out.println("Timer (phylogeny) started. (time: "
				+ (System.nanoTime() - starttime) * 1e-9 + " s)");
		}
		if (verbose) {
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out.println("Processing input...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: processing input");
		}
		List<int[]> columns = new ArrayList<int[]>();
		for (int i = 0; i < columns_char.size(); i++) {
			columns.add(MatrixTools.convertColumn(columns_char.get(i)));
		}
		if (verbose) {
			System.out.println("User wish for number of divisions: "
					+ userjobsnr);
		}

		// correct user input
		if (userjobsnr > length) {
			userjobsnr = length;
		} else if (userjobsnr < 1) {
			userjobsnr = 1;
		}
		int nrjobs = userjobsnr;

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		if (verbose) {
			System.out
					.println("Actual nr. of divisions in phylogenetic calculations: "
							+ nrjobs);
		}
		if (verbose) {
			System.out.println("Generating matrices for each node...");
		}

		tree.generateLeafList(names); // to speed up calculations later

		final double[][] probmatrix = new double[length][length];
		List<PhyloJob> jobs = new ArrayList<PhyloJob>();
		// First generate "exp(Rt)" matrix for each node.
		// This is the same for all columns of that node.
		tree.getRoot().calculateChildrenMatrix(param.getSD(), param.getSV(),
				param.getSV1());

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out.println("Generating jobs...");
		}
		if (verbose) {
			act.setCurrentActivity("Evolutionary model: dividing tasks");
		}

		// Create jobs for SINGLE columns
		// now all single columns are in one job
		// Add the only single-column job
		int colcnt = 0;
		PhyloJob lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob (for
		// debugging only)
		lastjob.startcol = 0;
		lastjob.jobid = 0;
		lastjob.endcol = length - 1;
		lastjob.type = false;
		lastjob.param = param;
		for (int col = 0; col < length; col++) {
			// send columns
			lastjob.columns.add(columns.get(col));
			colcnt++;
		}
		jobs.add(lastjob);

		// Create jobs for column PAIRS
		// count total number of column pairs
		int paircnt = 0;
		for (int i = 0; i < length; i++) {
			for (int j = i + 1; j < length; j++) {
				paircnt++;
			}
		}
		if(verbose){
			System.out.println("Total number of pairs: " + paircnt);
			System.out.println("Pairs in a job: " + (paircnt) / nrjobs);
		}
		// have to calculate matrices, but only once
		tree.getRoot().calculateChildrenMatrix(param.getDD(), param.getDV(),
				param.getDV1());

		int lastendcol = -1;
		int paircolcnt = 0;
		// add all jobs except last
		for (int jobnr = 0; jobnr < nrjobs - 1; jobnr++) {
			PhyloJob job = new PhyloJob();
			job.tree = Tree.copyTree(tree); // must copy the entire tree into
			// each phylojob
			job.names = names; // must copy the names into each phylojob
			job.type = true;
			job.param = param;
			job.jobid = jobnr + 1;

			int col1 = lastendcol + 1;
			job.startcol = col1;
			paircolcnt += length - col1;
			col1--;
			while (paircolcnt < ((long)(jobnr + 1) * (paircnt)) / nrjobs) {
				col1++;
				if (col1 == length) { // no more columns left, finish
					col1--;
					break;
				}
				job.columns.add(columns.get(col1));
				paircolcnt += length - col1;
			}

			job.endcol = col1;
			lastendcol = col1;
			for (int col2 = job.startcol + 1; col2 < length; col2++) {
				// send pairing columns
				job.columns2.add(columns.get(col2));
			}
			jobs.add(job);
			paircolcnt -= length - col1;
		}
		// Add the last job
		lastjob = new PhyloJob();
		lastjob.tree = Tree.copyTree(tree); // must copy the entire tree into
		// each phylojob
		lastjob.names = names; // must copy the names into each phylojob
		lastjob.startcol = lastendcol + 1;
		lastjob.endcol = length;
		lastjob.type = true;
		lastjob.param = param;
		lastjob.jobid = nrjobs;
		for (int col1 = lastjob.startcol; col1 < length; col1++) {
			// send columns
			lastjob.columns.add(columns.get(col1));
		}
		for (int col2 = lastjob.startcol + 1; col2 < length; col2++) {
			// send pairing columns
			lastjob.columns2.add(columns.get(col2));
		}
		jobs.add(lastjob);

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}

		if (verbose) {
			System.out
					.println("Total number of jobs in phylogenetic calculations: "
							+ jobs.size());
		}
		if (verbose) {
			System.out.println("Executing jobs...");
		}

		// start executing the jobs
		// counts how many phylojobs are done		
		final AtomicInteger finishedphylojobscount = new AtomicInteger(0); 

		act.setCurrentActivity("Evolutionary model: " +
					"calculating column probabilities");
		
		long gridstarttime = System.nanoTime();
		// execute single column jobs
		Progress singleColAct = act.getChildProgress(0.1);
		for (int jobnr = 0; jobnr < 1; jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;

			final Progress jobAct = singleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) colcnt);
			executor.startExecution(job, new JobListener() {

				public void jobFinished(double[][] result) {
					for (int col = 0; col < result.length; col++) {
						probmatrix[col + startcol][col + startcol] = result[col][0];
					}
					finishedphylojobscount.incrementAndGet();
					jobAct.setProgress(1.0);
				}
			});
		}
		
		/*for(int i = 0; i < probmatrix.length; i++){
			for(int j = 0; j < probmatrix[0].length; j++){
				System.out.print(probmatrix[i][j] + "\t");
			}
			System.out.println(); 
		}*/

		Progress doubleColAct = act.getChildProgress(0.9);

		// execute double column jobs
		for (int jobnr = 1; jobnr < jobs.size(); jobnr++) {
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
			PhyloJob job = jobs.get(jobnr);
			//final int jn = jobnr;
			final int startcol = job.startcol;
			final int endcol = job.endcol;
			final Progress jobAct = doubleColAct
					.getChildProgress((float) job.columns.size()
							/ (float) length);
			if (job.columns2.size() != 0) {
				executor.startExecution(job, new JobListener() {
					public void jobFinished(double[][] result) {
						for (int col1 = startcol; col1 < endcol + 1; col1++) {
							for (int col2 = col1 + 1; col2 < length; col2++) {
								probmatrix[col1][col2] = result[col1 - startcol][col2
										- startcol - 1];
								probmatrix[col2][col1] = probmatrix[col1][col2]; // make
								// it
								// symmetric
							}
						}
						finishedphylojobscount.incrementAndGet();
						jobAct.setProgress(1.0);
					}
				});
			} else {
				finishedphylojobscount.incrementAndGet();
				jobAct.setProgress(1.0);
			}
		}

		// wait for last thread to finish
		// wait for jobs to finish
		while (finishedphylojobscount.get() < jobs.size()) {
			Thread.sleep(100);
			if(act.shouldStop()){
				executor.shutDown();
			}
			act.checkStop();
		}
		singleColAct.setProgress(1.0);
		doubleColAct.setProgress(1.0);
		act.setProgress(1.0);
		//act.setCurrentActivity("Evolutionary model: done");

		if (verbose) {
			System.out.println("Done. (time: "
					+ (System.nanoTime() - starttime) * 1e-9 + " s)");
			System.out.println("Memory allocated: "
					+ (Runtime.getRuntime().totalMemory() / 1048576) + " MB");
			System.out.println("Memory used: "
					+ ((Runtime.getRuntime().totalMemory() - Runtime
							.getRuntime().freeMemory()) / 1048576) + " MB ");
		}
		
		if(verbose){


			System.out.println("TOTAL TIME ELAPSED IN PHYLOGENETIC PART: "
					+ (int)((System.nanoTime() - starttime) * 1e-9) + " seconds ");
			System.out.println("                    ...of which distributed: "
				+ (System.nanoTime() - gridstarttime) * 1e-9 + " seconds");
		}

		//System.out.println();
		

		// Result contains the a priori probability distribution matrix.
		
		return probmatrix;
	}
}