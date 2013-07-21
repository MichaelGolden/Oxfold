package uk.ac.ox.osscb.inoutside;



import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
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
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;
import uk.ac.ox.osscb.visualisation.DataVisualiser;

/**
 * Inside outside implementation
 * 
 * @author Vladimir
 */
public class CoFoldInsideOutsideCalculator {
	
	public static class ParallelInfo
	{
		private int upto = 0;
		List<Sector> sectors;
		
		public ParallelInfo(List<Sector> sectors)
		{
			this.sectors = sectors;
		}
		
		public synchronized Sector getMostCompleteSector()
		{
			for(int i = upto+1 ; i < sectors.size() ; i++)
			{
				if(sectors.get(i).finished)
				{
					upto = i;
				}
				else
				{
					break;
				}
			}
			
			return sectors.get(upto);
		}
		
		
		public synchronized boolean canRun(Sector s)
		{
			Sector upto = getMostCompleteSector();
			if(s.startb == upto.startb || (upto.startb+1 >= s.startb && s.endj  < upto.endj) || s.startb-1 == 0 || upto.n + 1 == s.n)
			{
				return true;
			}
			return false;
		}
	}
	
	public static class Sector
	{
		boolean finished;
		boolean running;
		int n;
		int startb;
		int endb;
		int startj;
		int endj;
		int starth;
		int endh;
		final InsideOutsideProbabilities iProbs;
		final InsideOutsideProbabilities oProbs;
		final NucleotideProbsPrecise pairingProbs;
		double alpha;
		double tau;
		int [] structure;
		boolean [][] canPair;
		//int [][] calculated;
		//int [][] calculates;
		ParallelInfo info;
		
		public Sector(int startb, int endb, int startj, int endj)
		{
			super();
			this.startb = startb;
			this.endb = endb;
			this.startj = startj;
			this.endj = endj;
			this.iProbs = null;
			this.oProbs = null;
			this.pairingProbs = null;
		}
		
		
		public Sector(int n, int startb, int endb, int startj, int endj, int starth, int endh, InsideOutsideProbabilities iProbs,
				 NucleotideProbsPrecise pairingProbs, double alpha, double tau,
				int[] structure, boolean[][] canPair, int [][] calculated, ParallelInfo info) {
			super();
			this.n = n;
			this.startb = startb;
			this.endb = endb;
			this.startj = startj;
			this.endj = endj;
			this.starth = starth;
			this.endh = endh;
			this.iProbs = iProbs;
			this.pairingProbs = pairingProbs;
			this.alpha = alpha;
			this.tau = tau;
			this.structure = structure;
			this.canPair = canPair;
			//this.calculated =calculated;
			this.oProbs = null;
			this.info = info;
			//this.calculates = new int[calculated.length][calculated.length];
		}
		
		public Sector(int n, int startb, int endb, int startj, int endj, int starth, int endh, InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
				NucleotideProbsPrecise pairingProbs, double alpha, double tau,
				int[] structure, boolean[][] canPair, int [][] calculated, ParallelInfo info) {
			super();
			this.n = n;
			this.startb = startb;
			this.endb = endb;
			this.startj = startj;
			this.endj = endj;
			this.starth = starth;
			this.endh = endh;
			this.iProbs = iProbs;
			this.oProbs = oProbs;
			this.pairingProbs = pairingProbs;
			this.alpha = alpha;
			this.tau = tau;
			this.structure = structure;
			this.canPair = canPair;
			//this.calculated =calculated;
			this.info = info;
			//this.calculates = new int[calculated.length][calculated.length];
		}



		@Override
		public String toString() {
			return "Sector [startb=" + startb + ", endb=" + endb + ", startj="
					+ startj + ", endj=" + endj+"]";
		}
	}
	
	private Grammar grammar;

	public CoFoldInsideOutsideCalculator(Grammar grammar) {
		super();

		if(null == grammar)
			throw new NullPointerException("grammar cannot be null");

		this.grammar = grammar;
	}
	
	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#inside(uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, int[], boolean[][])
	 */
	public InsideOutsideProbabilities inside(NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}	
		}
		
		return iProbs;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#inside(uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, int[], boolean[][])
	 */
	public InsideOutsideProbabilities insideWithVis(NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair) {

		
		
		final int seqLen = structure.length;
		
		double [][] calculation = new double[seqLen][seqLen];
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
					calculation[charIdx][charIdx] = 1.0;
				}
			}

			/*
			try {
				iProbs.writeTable(new File("iprobs_"+rule1.getLeft()+".txt"), rule1.getLeft());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/
			
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				for(int x = 0 ;  x < calculation.length ; x++)
				{
					for(int y = 0 ;  y < calculation.length ; y++)
					{
						if(calculation[x][y] != 0)
						{
							calculation[x][y]=0.4;
						}
					}
				}
				
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
						calculation[j][h] = 0.7;
						calculation[h+1][j+b] = 0.7;
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);
					
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
						calculation[j][j+b] = 0.7;
						calculation[j+1][j+b-1] = 0.7;
					}
				}
				
				calculation[j][j+b] = 1;
				
				try {
					DecimalFormat df = new DecimalFormat("000");
					DataVisualiser.drawBasePairProbMatrix(calculation, 400, 400).savePNG(new File("algorithm/in_"+df.format(b)+"_"+df.format(j)+".svg"), new File("algorithm/in_"+df.format(b)+"_"+df.format(j)+".png"));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}	
			
		
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}
	
	public InsideOutsideProbabilities insideParallel ( NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair)
	{
		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}			
		}
		
		

		runInsideJobs(iProbs,pairingProbs,alpha,tau,structure,canPair);
		return iProbs;
	}
	
	public void runInsideJobs(InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{

		boolean run_parallel = true;
		/*int [][] calculated = new int[structure.length][structure.length];
		for(int i = 0 ; i < calculated.length ; i++)
		{
			calculated[i][i]++;
		}*/
		if(!run_parallel)
		{

		
			/*for(int b = 1 ; b < structure.length ; b++)
			{
				ArrayList<Sector> sectors2 = createInsideSectors(b, structure.length, 8, iProbs, pairingProbs,alpha,tau,structure,canPair, null);
				for(Sector s: sectors2)
				{
					runInsideJob(s);
				}
			}*/
		
		}
		else
		{
			
			try {
				
				//processInputs(sectors);
				/*for(int b = 1 ; b < structure.length ; b++)
				{
					ArrayList<Sector> sectors2 = createInsideSectors(b, structure.length, 8, iProbs, pairingProbs,alpha,tau,structure,canPair, calculated);
					processInsideSectors(sectors2);
				}*/
				ArrayList<Sector> sectors = createInsideSectors2(structure.length, iProbs, pairingProbs,alpha,tau,structure,canPair);
				processInsideSectors(sectors);
				/*for(int b = 1 ; b < structure.length ; b++)
				{
					ArrayList<Sector> sectors2 = createInsideSectors(b, structure.length, 8, iProbs, pairingProbs,alpha,tau,structure,canPair, calculated);
					
				}
				processInsideSectors(sectors2);*/
			
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public void runOutsideJobs(InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{

		boolean run_parallel = true;
		
		if(!run_parallel)
		{
			//ArrayList<Sector> sectors2 = createOutsideSectors(0, structure.length, 8, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair);
			//runOutsideJob(sectors2.get(0), sectors2);
			
			/*
			for (int b = structure.length-2; b >= 0; b--) 
			{
				//ArrayList<Sector> sectors2 = createOutsideSectors(b, structure.length, 8, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair);
				for(Sector s: sectors2)
				{
					runOutsideJob(s, sectors2);
					//break;
				}
			}*/
		}
		else
		{
			
			try {
				
				//System.out.println("Start outside creation");
				//long t1 = System.nanoTime();
				ArrayList<Sector> sectors = createOutsideSectors2(structure.length, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair);
				//long t2 = System.nanoTime();
				//System.out.println(((t2-t1)/1e9)+"s");
				//System.out.println("End outside creation "+sectors.size());
				processOutsideSectors(sectors);
				/*for (int b = structure.length-2; b >= 0; b--) 
				{
					ArrayList<Sector> sectors2 = createOutsideSectors(b, structure.length, 8, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair);
					processOutsideSectors(sectors2);
				}*/
			
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	
	public ArrayList<Sector> createInsideSectors2(int seqLen, InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{
		int divisions = threads*2;
		ArrayList<Sector> sectors = new ArrayList<Sector>();
		ParallelInfo info = new ParallelInfo(sectors);
		for (int b = 1; b <= seqLen; b += 1) {
			int j_divisions = Math.max(2, (seqLen-b)/divisions);			
			for (int j = 0; j <= seqLen-b; j += j_divisions) {
				Sector s = new Sector(sectors.size(), b, Math.min(seqLen, b+1), j,Math.min(seqLen-b, j+j_divisions), 0, seqLen, iProbs, pairingProbs,alpha,tau,structure,canPair, null, info);
				sectors.add(s);
			}
		}
		
		return sectors;
	}
	
	public ArrayList<Sector> createOutsideSectors2(int seqLen, InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{
		int divisions = threads*2;
		ArrayList<Sector> sectors = new ArrayList<Sector>();
	

		ParallelInfo info = new ParallelInfo(sectors);
		for (int b = structure.length-2; b >= 0; b--) 
		{
			int j_divisions = Math.max(2, (seqLen-b)/divisions);
			for (int j = 0; j <= seqLen-b; j += j_divisions) {
				Sector s = new Sector(sectors.size(), b, Math.max(0, b-1), j,Math.min(seqLen-b, j+j_divisions), 0, seqLen, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair, null, info);
				sectors.add(s);
			}
		}
		
		return sectors;
	}
	
	public void runInsideJob(Sector sector) {

    	ParallelInfo info = sector.info;
    	while(!info.canRun(sector))
    	{
    		
    	}
		
		InsideOutsideProbabilities iProbs = sector.iProbs;
		NucleotideProbsPrecise pairingProbs = sector.pairingProbs;
		boolean [][] canPair = sector.canPair;
		double alpha = sector.alpha;
		double tau = sector.tau;
		for (int b = sector.startb; b < sector.endb ; b++) {
			for (int j = sector.startj ; j < sector.endj ; j++) {
			//System.out.println(b+"\t"+j+"\t");
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j ; h < j+b  ; h++) {
						/*if(sector.calculated[j][h] < 1)
						{
							System.out.println("AWARNING "+j+","+h+" not calculated "+sector);					
							int j1 = j;
							int b1 = h-j1;
							System.out.println("Relies "+b1+"\t"+j1);
							System.exit(0);
						}
						if(sector.calculated[h+1][j+b] < 1)
						{
							System.out.println("BWARNING diag="+b+"\tpos="+j+"\t"+(h+1)+","+(j+b)+" not calculated "+sector);
							int j1 = h+1;
							int b1 = (j+b)-j1;
							System.out.println("Relies "+b1+"\t"+j1);
							System.exit(0);
						}*/
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);

					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);
					

					//sector.calculated[j][j+b]++;
				}
				
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}

				//sector.calculated[j][j+b]++;
			}	
			
		
		}
	}
	

	
	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	public InsideOutsideProbabilities outsideParallel(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;		
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		runOutsideJobs(insideProbs, oProbs,pairingProbs,alpha,tau,structure,canPair);
	
		
		return oProbs;
	}
	

	
	public void runOutsideJob(Sector sector, List<Sector> sectors) {
		
		ParallelInfo info = sector.info;
    	while(!info.canRun(sector))
    	{
    		
    	}
		int seqLen = sector.structure.length;
		
		InsideOutsideProbabilities insideProbs = sector.iProbs;
		InsideOutsideProbabilities oProbs = sector.oProbs;
		NucleotideProbsPrecise pairingProbs = sector.pairingProbs;
		boolean [][] canPair = sector.canPair;
		double alpha = sector.alpha;
		double tau = sector.tau;
		
		int b = sector.startb;
		// iterate over (decreasing) length of inside interval
		//for (int b = seqLen-2; b >= 0; b--) {
		// iterate over (decreasing) length of inside interval
			//for (int b = sector.startb ; b >= sector.endb; b--) {
				for (int j = sector.startj; j < sector.endj ; j++) {
					//for (int j = 0; j < seqLen - b; j++) {
						for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
							PointRes tmp = PointRes.ZERO;
							for (int k = j + b + 1; k < seqLen; k++) {
								PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
								PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
								tmp = tmp.add(oProb.multiply( iProb));
							}
							PointRes probIncrement = tmp.multiply(rule3.getProbability());
							oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
							tmp = PointRes.ZERO;
							for(int k = 0; k<j;k++) {
								PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
								PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
								tmp = tmp.add(oProb.multiply(iProb));						
							}
							probIncrement = tmp.multiply(rule3.getProbability());					
							oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
						}
						// get contributions of rules of type U->(V)
						if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
							for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
								PointRes probIncrement = PointRes.ZERO;
								double distance = Math.abs(b);
								PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
								//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
								probIncrement = rule2.getProbability().multiply(exp)
										.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
										.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
								oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
							}
						}
					}
			//}
	}
	


    final int threads = 4;
	 public List<Object> processInsideSectors(final List<Sector> inputs)
	            throws InterruptedException, ExecutionException {

	        ExecutorService service = Executors.newFixedThreadPool(threads);

	        List<Future<Object>> futures = new ArrayList<Future<Object>>();
	        for (final Sector input : inputs) {
	            Callable<Object> callable = new Callable<Object>() {
	                public Object call() throws Exception {
	                	input.running  = true;
                		//runInsideJob(input);
	                	runInsideJob(input);
                		//System.out.println("Current:"+input+"\t MR"+input.info.getMostCompleteSector());
                		input.finished = true;
	                    return "Finished";
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
	 
	 public List<Object> processInsideSectors3(final List<Sector> inputs)
	            throws InterruptedException, ExecutionException {

	        final int threads = 4;
	        ExecutorService service = Executors.newFixedThreadPool(threads);

	        List<Future<Object>> futures = new ArrayList<Future<Object>>();
	        for (final Sector input : inputs) {
	            Callable<Object> callable = new Callable<Object>() {
	                public Object call() throws Exception {
	                	input.running  = true;
             		runInsideJob(input);  
             		input.finished = true;
	                    return "Finished";
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
	 
	 public List<Object> processOutsideSectors(final List<Sector> inputs)
	            throws InterruptedException, ExecutionException {

	        ExecutorService service = Executors.newFixedThreadPool(threads);

	        List<Future<Object>> futures = new ArrayList<Future<Object>>();
	        for (final Sector input : inputs) {
	            Callable<Object> callable = new Callable<Object>() {
	                public Object call() throws Exception {
	                	input.running  = true;
             		runOutsideJob(input, inputs);  
             		input.finished = true;
	                    return "Finished";
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


	public InsideOutsideProbabilities insideE(NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}			
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}
	
	
	
	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	public InsideOutsideProbabilities outside3(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;		
		double [][] calculation = new double[seqLen][seqLen];
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		calculation[0][seqLen-1] = 1;
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {

				for(int x = 0 ;  x < calculation.length ; x++)
				{
					for(int y = 0 ;  y < calculation.length ; y++)
					{
						if(calculation[x][y] != 0)
						{
							calculation[x][y]=0.4;
						}
					}
				}
				
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
						
						calculation[j][k]= 0.7;
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));	
						
						calculation[k][j+b]= 0.7;
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						double distance = Math.abs(b);
						PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
						//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability().multiply(exp)
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
						
						calculation[j-1][j+b+1]= 0.7;
					}
				}
				

				
				calculation[j][j+b]=1;
				
				boolean vis = false;
				if(vis)
				{
				try {
					DecimalFormat df = new DecimalFormat("000");
					System.out.println("here");
					int bb = seqLen+2-b;
					DataVisualiser.drawBasePairProbMatrix(calculation, 400, 400).savePNG(new File("algorithm/out_"+df.format(bb)+"_"+df.format(j)+".svg"), new File("algorithm/out"+df.format(bb)+"_"+df.format(j)+".png"));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} 
				}
			}
		}
		
		return oProbs;
	}
	
	public InsideOutsideProbabilities outsideE(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
	}

	//@Override
	public InsideOutsideProbabilities inside(
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double alpha, double tau, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		
		InsideOutsideProbabilities iProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise inside with values of i[U,j,j]
		for (ProductionRule rule1 : grammar.getRules(RuleType.RULE1)) {
			for (int charIdx = 0; charIdx<seqLen; charIdx++) {
				if(Constants.UnpairedBaseIdx == structure[charIdx]) {
					PointRes prob = pairingProbs.getUnpairingProbability(charIdx).multiply(rule1.getProbability());
					iProbs.setProb(rule1.getLeft(), charIdx, charIdx, prob);
				}
			}
		}
		
		// iterate over length of inside interval
		for (int b = 1; b < seqLen; b++) {
			for (int j = 0; j < seqLen-b; j++) {
				// deal with contribution of rules of type U->VW
				for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int h = j; h < j+b; h++) {
						PointRes prob1 = iProbs.getProb(rule3.getRight()[0], j, h);
						PointRes prob2 = iProbs.getProb(rule3.getRight()[1], h+1, j+b);
						tmp = tmp.add(prob1.multiply(prob2));
					}
					PointRes probIncrement = rule3.getProbability().multiply(tmp);
					iProbs.increment(rule3.getLeft(), j, j+b, probIncrement);						
				}
				// deal with production rules of type U->(V)				
				if (canPair[j][j+b]) {
					for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
						double distance = Math.abs(b);
						double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
						PointRes probIncrement = PointRes.ZERO;
						probIncrement = rule2.getProbability().multiply(PointRes.valueOf(exp))
								.multiply(pairingProbs.getPairingProbability(j, j+b))
								.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
						iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
					}
				}
			}			
		}
		
		// System.out.println(iProbs.printAllTables());
		
		return iProbs;
	}

	//@Override
	public InsideOutsideProbabilities outside(
			InsideOutsideProbabilities insideProbs,
			NucleotideProbsPrecise pairingProbs, double[][] distances,
			double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
		
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						double distance = Math.abs(b);
						PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
						probIncrement = rule2.getProbability().multiply(exp)
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
		
	}
	
	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	public InsideOutsideProbabilities outside(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair) {
		
		final int seqLen = structure.length;
				
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		
		// initialise outside probabilities
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		
		// iterate over (decreasing) length of inside interval
		for (int b = seqLen-2; b >= 0; b--) {
			for (int j = 0; j < seqLen - b; j++) {
				for (ProductionRule rule3 : this.grammar.getRules(RuleType.RULE3)) {
					PointRes tmp = PointRes.ZERO;
					for (int k = j + b + 1; k < seqLen; k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), j, k);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[1], j+b+1, k);
						tmp = tmp.add(oProb.multiply( iProb));
					}
					PointRes probIncrement = tmp.multiply(rule3.getProbability());
					oProbs.increment(rule3.getRight()[0], j, j+b, probIncrement);
					tmp = PointRes.ZERO;
					for(int k = 0; k<j;k++) {
						PointRes oProb = oProbs.getProb(rule3.getLeft(), k, j+b);
						PointRes iProb = insideProbs.getProb(rule3.getRight()[0], k, j-1);
						tmp = tmp.add(oProb.multiply(iProb));						
					}
					probIncrement = tmp.multiply(rule3.getProbability());					
					oProbs.increment(rule3.getRight()[1], j, j+b, probIncrement);
				}
				// get contributions of rules of type U->(V)
				if ((j>=1)&&(j+b+1<seqLen)&&(canPair[j-1][j+b+1])) {
					for (ProductionRule rule2: this.grammar.getRules(RuleType.RULE2)) {
						PointRes probIncrement = PointRes.ZERO;
						double distance = Math.abs(b);
						PointRes exp = PointRes.valueOf(alpha*(Math.exp(-distance/tau) - 1) + 1);
						//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability().multiply(exp)
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		
		return oProbs;
	}
}
