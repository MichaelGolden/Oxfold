package uk.ac.ox.osscb.inoutside;


import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;

public class PosteriorProbabilitiesCalculator {
		
		private Grammar grammar;
		
		public PosteriorProbabilitiesCalculator(Grammar grammar) {
			this.grammar = grammar;
		}
		
		public PosteriorProbabilities calculate(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				int[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							pairedProbs[j][k] = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);

							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}


		public ShortDoublePosteriorProbabilities calculateShort(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			double[][] pairedProbs = new double[length][length];
			double[] unpairedProbs = new double[length];
			
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					PointRes prob = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
					unpairedProbs[j] = prob.doubleValue();
				}
			}
			//paired probabilities
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							PointRes prob = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
							pairedProbs[j][k] = prob.doubleValue();
						}
					}
				}
			}
			ShortDoublePosteriorProbabilities posteriorProbabilities = new ShortDoublePosteriorProbabilities(unpairedProbs, pairedProbs);
			return posteriorProbabilities;
		}
		
		public PosteriorProbabilities calculateE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							pairedProbs[j][k] = tmp.multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
													
							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}
		
		public static class PosteriorProbabilitiesJob
		{
			InsideOutsideProbabilities insideProb;
			InsideOutsideProbabilities outsideProbs;
			NucleotideProbsPrecise nucleotideProbs;
			int istart;
			int iend;
			double[][] distances;
			double weight;
			int[] structure;
			boolean[][] canPair;
			PointRes[][] pairedProbs;
			PointRes[] unpairedProbs;
			public PosteriorProbabilitiesJob(int istart, int iend,
					InsideOutsideProbabilities insideProb,
					InsideOutsideProbabilities outsideProbs,
					NucleotideProbsPrecise nucleotideProbs,
					double[][] distances, double weight, int[] structure,
					boolean[][] canPair, PointRes[][] pairedProbs, PointRes[] unpairedProbs) {
				super();
				this.istart = istart;
				this.iend = iend;
				this.insideProb = insideProb;
				this.outsideProbs = outsideProbs;
				this.nucleotideProbs = nucleotideProbs;
				this.distances = distances;
				this.weight = weight;
				this.structure = structure;
				this.canPair = canPair;
				this.pairedProbs = pairedProbs;
				this.unpairedProbs = unpairedProbs;
			}
			@Override
			public String toString() {
				return "PosteriorProbabilitiesJob [istart=" + istart
						+ ", iend=" + iend + "]";
			}		
			
			
		}
		
		public  static class PosteriorProbabilitiesResult
		{
			PointRes maxp = PointRes.ZERO; 
			int rightIdx = -1;
			int leftIdx = -1;
		}
		
		public ArrayList<PosteriorProbabilitiesJob> createJobs(int threads, InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair)
			{
				ArrayList<PosteriorProbabilitiesJob> jobs = new ArrayList<PosteriorProbabilitiesJob>();

				PointRes [][] pairedProbs = new PointRes[structure.length][structure.length];
				PointRes [] unpairedProbs = new PointRes[structure.length];
				int divisions = Math.max(2, structure.length/threads/2);
				for(int i = 0 ; i <= structure.length ; i += divisions)
				{
					jobs.add(new PosteriorProbabilitiesJob(i, Math.min(structure.length, i+divisions), insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair, pairedProbs, unpairedProbs));
				}
				//System.out.println(jobs+"\t"+structure.length);
				return jobs;
			}
		

	    public static final int threads = Constants.threads;;
		public List<PosteriorProbabilitiesResult> processInputs(List<PosteriorProbabilitiesJob> inputs)
		        throws InterruptedException, ExecutionException {

		    ExecutorService service = Executors.newFixedThreadPool(threads);

		    List<Future<PosteriorProbabilitiesResult>> futures = new ArrayList<Future<PosteriorProbabilitiesResult>>();
		    for (final PosteriorProbabilitiesJob input : inputs) {
		        Callable<PosteriorProbabilitiesResult> callable = new Callable<PosteriorProbabilitiesResult>() {
		            public PosteriorProbabilitiesResult call() throws Exception {
		            	PosteriorProbabilitiesResult output = calculate(input);
		                return output;
		            }
		        };
		        futures.add(service.submit(callable));
		    }

		    service.shutdown();

		    List<PosteriorProbabilitiesResult> outputs = new ArrayList<PosteriorProbabilitiesResult>();
		    for (Future<PosteriorProbabilitiesResult> future : futures) {
		        outputs.add(future.get());
		    }
		    return outputs;
		}
		
		public PosteriorProbabilities calculateParallel(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair) {
			ArrayList<PosteriorProbabilitiesJob> jobs = createJobs(threads, insideProbs, outsideProbs, nucleotideProbs, distances, weight, structure, canPair);
			try {
				List<PosteriorProbabilitiesResult> results = processInputs(jobs);
				
				PointRes maxp = null;
				int leftIdx = -1;
				int rightIdx = -1;
				for(PosteriorProbabilitiesResult result : results)
				{
					if (maxp == null || result.maxp.compareTo(maxp) > 0) {
						maxp = result.maxp; 
						leftIdx = result.leftIdx;
						rightIdx = result.rightIdx;
					}
				}
				//System.out.println("K"+maxp+"\t"+leftIdx+"\t"+rightIdx);
				
				PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(jobs.get(0).unpairedProbs, jobs.get(0).pairedProbs, maxp, leftIdx, rightIdx);
				/*for(int i = 0 ; i < posteriorProbabilities.pairedProbs.length ; i++)
				{
					
					if(posteriorProbabilities.unpairedProbs[i] == null)
					{
						System.err.println(i+"\t"+posteriorProbabilities.unpairedProbs[i]+"\t");
					}
					for(int j = 0 ; j < posteriorProbabilities.pairedProbs.length ; j++)
					{
						if(posteriorProbabilities.pairedProbs[i][j] == null)
						{
							System.err.println(i+"\t"+j+"\t"+posteriorProbabilities.pairedProbs[i][j]+"\t");
						}
					}
				}*/
				return posteriorProbabilities;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			return null;
		}
		
		public PosteriorProbabilitiesResult calculate(PosteriorProbabilitiesJob job)
		{
			int [] structure = job.structure;
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = job.pairedProbs;
			PointRes[] unpairedProbs = job.unpairedProbs;
			InsideOutsideProbabilities insideProbs = job.insideProb;
			InsideOutsideProbabilities outsideProbs = job.outsideProbs;
			NucleotideProbsPrecise nucleotideProbs = job.nucleotideProbs;
			boolean [][] canPair = job.canPair;
			double [][] distances =job.distances;
			double weight = job.weight;
			//just zeroing the arrays

			for(int i = job.istart ; i < job.iend ; i++)
			{
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			
			if(job.istart == 0)
			{
				//unpaired probabilities	
				for (int j=0; j<length; j++) {

					unpairedProbs[j] = PointRes.ZERO;
					if (Constants.UnpairedBaseIdx == structure[j]) {
						PointRes tmp = PointRes.ZERO;
						for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
							tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
						}
						unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
								.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
					}
				}
			}
			
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=job.istart; j<job.iend; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							pairedProbs[j][k] = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);

							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			
			PosteriorProbabilitiesResult res = new PosteriorProbabilitiesResult();
			res.maxp = maxp;
			res.leftIdx = leftIdx;
			res.rightIdx = rightIdx;
			return res;
		}
		
		//takes double[][] for distance param (not int[][])
		public PosteriorProbabilities calculate2(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs,
				double[][] distances, double weight, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			PointRes maxp = PointRes.ZERO; 
			PointRes[][] pairedProbs = new PointRes[length][length];
			PointRes[] unpairedProbs = new PointRes[length];
			//just zeroing the arrays
			for (int i=0; i<length; i++) {
				unpairedProbs[i] = PointRes.ZERO;
				for (int j=0; j<length; j++) {
				pairedProbs[i][j] = PointRes.ZERO;	
				}
			}
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					unpairedProbs[j] = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
				}
			}
			//paired probabilities
			int rightIdx = -1;
			int leftIdx = -1;
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes tmpExp = PointRes.valueOf(Math.exp(-distances[j][k]/weight));
							pairedProbs[j][k] = tmp.multiply(tmpExp).multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);

							if (pairedProbs[j][k].compareTo(maxp) > 0) {
								maxp = pairedProbs[j][k]; 
								leftIdx = j;
								rightIdx = k;
							}
						}
					}
				}
			}
			PosteriorProbabilities posteriorProbabilities = new PosteriorProbabilities(unpairedProbs, pairedProbs, maxp, leftIdx, rightIdx);
			return posteriorProbabilities;
		}
		
		public ShortDoublePosteriorProbabilities calculateShortE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, int[] structure, boolean[][] canPair) {
			int length = structure.length;
			double[][] pairedProbs = new double[length][length];
			double[] unpairedProbs = new double[length];
			
			//unpaired probabilities	
			for (int j=0; j<length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					PointRes tmp = PointRes.ZERO;
					for(ProductionRule pr : this.grammar.getRules(RuleType.RULE1)){
						tmp = tmp.add(outsideProbs.getProb(new Character(pr.getLeft()), j, j).multiply(pr.getProbability()));
					}
					PointRes prob = tmp.multiply(nucleotideProbs.getUnpairingProbability(j))
							.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
					unpairedProbs[j] = prob.doubleValue();
				}
			}
			//paired probabilities
			for(int j=0; j<length; j++){
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for(int k=j+4; k<length; k++){
						if (canPair[j][k]) {
							PointRes tmp = PointRes.ZERO;
							for(ProductionRule pr : this.grammar.getRules(RuleType.RULE2)){
								PointRes inside = insideProbs.getProb(new Character(pr.getRight()[1]), j+1, k-1);
								PointRes outside = outsideProbs.getProb(new Character(pr.getLeft()), j, k);
								PointRes addition = outside.multiply(pr.getProbability()).multiply(inside);
								tmp = tmp.add(addition);
							}
							PointRes prob = tmp.multiply(nucleotideProbs.getPairingProbability(j, k))
									.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, length-1), RoundingMode.HALF_UP);
							pairedProbs[j][k] = prob.doubleValue();
						}
					}
				}
			}
			ShortDoublePosteriorProbabilities posteriorProbabilities = new ShortDoublePosteriorProbabilities(unpairedProbs, pairedProbs);
			return posteriorProbabilities;
		}
		
		public PointRes calculateCompE(boolean[][] incomp, Helix helix, InsideOutsideProbabilities insideProbs, InsideOutsideProbabilities outsideProbs,
				NucleotideProbsPrecise nucleotideProbs, PointRes[][] pairedProbs) {
			PointRes rprob = PointRes.ZERO;
			for (int j=0; j<incomp.length; j++) {
				for (int k = j+1; k<incomp.length; k++) {
					if (incomp[j][k]) {
						rprob = rprob.add(pairedProbs[j][k]);
					}
				}
			}
			PointRes hprob = findHelixProbabilityE(insideProbs, outsideProbs, nucleotideProbs, helix);
			if (rprob.signum()>0) {
				rprob = hprob.divide(rprob.add(hprob), RoundingMode.HALF_UP);
			}
			return rprob;
		}
		
		public static class Interval
		{
			int start;
			int end;
			PointRes[][] pairedProbs;
			PointRes[] unpairedProbs;
			boolean[][] canPair;
			int length;
			PointRes [][] diffs;
			public Interval(int start, int end, PointRes[][] pairedProbs,
					PointRes[] unpairedProbs, boolean[][] canPair, int length,
					PointRes[][] diffs) {
				super();
				this.start = start;
				this.end = end;
				this.pairedProbs = pairedProbs;
				this.unpairedProbs = unpairedProbs;
				this.canPair = canPair;
				this.length = length;
				this.diffs = diffs;
			}	
			
			
		}
		
		public List<Object> processInputsForDiffs(List<Interval> inputs)
		        throws InterruptedException, ExecutionException {

		    ExecutorService service = Executors.newFixedThreadPool(threads);

		    List<Future<Object>> futures = new ArrayList<Future<Object>>();
		    for (final Interval input : inputs) {
		        Callable<Object> callable = new Callable<Object>() {
		            public Object call() throws Exception {
		            	
		            	calculateDiffs(input);
		                return null;
		            }
		        };
		        futures.add(service.submit(callable));
		    }

		    service.shutdown();

		    List<Object> outputs = new ArrayList<Object>();
		    for (Future<Object> future : futures) {
		        outputs.add(future.get());
		    }
		    return outputs;
		}
		
		public void calculateDiffs(Interval interval)
		{
			PointRes [][] pairedProbs = interval.pairedProbs;
			PointRes [] unpairedProbs = interval.unpairedProbs;
			boolean [][] canPair = interval.canPair;
			int length = interval.length;
			PointRes [][] diffs = interval.diffs;			
			
			for(int j = interval.start ; j < interval.end ; j++)
        	{
            	for (int k = 0; k<length; k++) {
					if (canPair[j][k]) {
						diffs[j][k] = pairedProbs[j][k].subtract((unpairedProbs[j].add(unpairedProbs[k])).divide(PointRes.valueOf(2)));
					} else {
						diffs[j][k] = PointRes.valueOf(-1);
					}	
				}
        	}
		}
		
		/**
		 * calculate accuracy scores
		 */
		public PointRes[][] getDiffs(final PointRes[][] pairedProbs, final PointRes[] unpairedProbs, final boolean[][] canPair) {
			/*int length = unpairedProbs.length;
			PointRes[][] diffs = new PointRes[length][length];
			for (int j = 0; j<length; j++) {
				for (int k = 0; k<length; k++) {
					if (canPair[j][k]) {
						diffs[j][k] = pairedProbs[j][k].subtract((unpairedProbs[j].add(unpairedProbs[k])).divide(PointRes.valueOf(2)));
					} else {
						diffs[j][k] = PointRes.valueOf(-1);
					}	
				}
			}
			return diffs;
			*/
			final int length = unpairedProbs.length;
			final PointRes[][] diffs = new PointRes[length][length];
			ArrayList<Interval> intervals = new ArrayList<Interval>();
			int divisions = Math.max(2, length/threads/2);
			for(int i = 0 ; i <= length ; i += divisions)
			{
				intervals.add(new Interval(i, Math.min(length, i+divisions), pairedProbs, unpairedProbs, canPair, length, diffs));
			}
			try {
				processInputsForDiffs(intervals);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return intervals.get(0).diffs;
		}
		
		
		/**
		 * calculate probability of helix forming
		 */
		
		private PointRes findHelixProbabilityE(InsideOutsideProbabilities insideProbs,
				InsideOutsideProbabilities outsideProbs, NucleotideProbsPrecise nucleotideProbs, Helix helix) {
			int left = helix.getLeftIdx(); int right = helix.getRightIdx();
			int length = helix.getHelixLength();
			PointRes nprob = PointRes.ONE;
			for (int j = 0; j<length; j++) {
				nprob = nprob.multiply(nucleotideProbs.getPairingProbability(left+j, right-j));
			}
			nprob = nprob.divide(insideProbs.getProb(this.grammar.getNonterminals()[0], 0, insideProbs.getDimension()-1), RoundingMode.HALF_UP);
			HashMap<Character,PointRes> Probs = new HashMap<Character,PointRes>();
			//initialise Hashmap with inside probabilities
			for (ProductionRule pr : this.grammar.getRules(RuleType.RULE2)) {
				Probs.put(pr.getLeft(), insideProbs.getProb(pr.getRight()[1], left+length, right-length).multiply(pr.getProbability())); 
			}
			for (int j = length-2; j>=0; j--) {
				HashMap<Character,PointRes> tmpProbs = new HashMap<Character,PointRes>();
				for (char nt: this.grammar.getNonterminals()) {
					tmpProbs.put(nt, PointRes.ZERO);
				}
				for (ProductionRule pr: this.grammar.getRules(RuleType.RULE2)) {
					PointRes newtmp = pr.getProbability().multiply(Probs.get(pr.getRight()[1]));
					tmpProbs.put(pr.getLeft(), newtmp);
				}
				Probs = tmpProbs;
			}
			PointRes prob = PointRes.ZERO;
			for (char nt : this.grammar.getNonterminals()) {
				prob = prob.add(outsideProbs.getProb(nt, left, right).multiply(Probs.get(nt)));
			}
			return prob;
		}
}
