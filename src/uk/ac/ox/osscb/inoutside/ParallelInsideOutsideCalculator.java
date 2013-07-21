package uk.ac.ox.osscb.inoutside;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
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

public class ParallelInsideOutsideCalculator {

	final int threads = Constants.threads;


	private Grammar grammar;

	public ParallelInsideOutsideCalculator(Grammar grammar) {
		super();

		if(null == grammar)
			throw new NullPointerException("grammar cannot be null");

		this.grammar = grammar;
	}
	


	public InsideOutsideProbabilities insideParallel (NucleotideProbsPrecise pairingProbs, int [] structure, boolean [][] canPair, PointResUpperMatrix basepairWeighting)
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
		runInsideJobs(iProbs,pairingProbs,structure,canPair, basepairWeighting);
		return iProbs;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	public InsideOutsideProbabilities outsideParallel(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting) {

		final int seqLen = structure.length;
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		runOutsideJobs(insideProbs, oProbs,pairingProbs,structure,canPair, basepairWeighting);

		return oProbs;
	}

	public InsideOutsideProbabilities insideParallelCoFold (NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair)
	{

		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				double distance = Math.abs(i-j);
				double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
				basepairWeighting.set(i, j, PointRes.valueOf(exp));
			}
		}

		return insideParallel (pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public InsideOutsideProbabilities outsideParallelCoFold (InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair)
	{

		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				double distance = Math.abs(i-j);
				double exp = alpha*(Math.exp(-distance/tau) - 1) + 1;
				basepairWeighting.set(i, j, PointRes.valueOf(exp));
			}
		}

		return outsideParallel (insideProbs, pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public InsideOutsideProbabilities insideParallelOxfold (NucleotideProbsPrecise pairingProbs,  int[][] distances, double weight, int [] structure, boolean [][] canPair)
	{
		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				basepairWeighting.set(i, j, PointRes.valueOf(Math.exp(-distances[i][j]/weight)));
			}
		}

		return insideParallel (pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public InsideOutsideProbabilities outsideParallelOxfold (InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, int[][] distances, double weight, int[] structure, boolean[][] canPair)
	{

		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				basepairWeighting.set(i, j, PointRes.valueOf(Math.exp(-distances[i][j]/weight)));
			}
		}
		return outsideParallel (insideProbs, pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public InsideOutsideProbabilities insideParallelOxfold (NucleotideProbsPrecise pairingProbs,  double[][] distances, double weight, int [] structure, boolean [][] canPair)
	{
		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				basepairWeighting.set(i, j, PointRes.valueOf(Math.exp(-distances[i][j]/weight)));
			}
		}

		return insideParallel (pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public InsideOutsideProbabilities outsideParallelOxfold (InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double[][] distances, double weight, int[] structure, boolean[][] canPair)
	{

		PointResUpperMatrix basepairWeighting = new PointResUpperMatrix(structure.length);
		for(int i = 0 ; i < structure.length ; i++)
		{
			for(int j = i ; j < structure.length ; j++)
			{
				basepairWeighting.set(i, j, PointRes.valueOf(Math.exp(-distances[i][j]/weight)));
			}
		}
		return outsideParallel (insideProbs, pairingProbs, structure,canPair, basepairWeighting);
	}
	
	public void runInsideJobs(InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs,
			int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting)
	{

		try {
			ArrayList<Sector> sectors = createInsideSectors(structure.length, iProbs, pairingProbs,structure,canPair, basepairWeighting);
			processInsideSectors(sectors);

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void runOutsideJobs(InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
			NucleotideProbsPrecise pairingProbs,int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting)
	{

		try {
			ArrayList<Sector> sectors = createOutsideSectors(structure.length, iProbs, oProbs, pairingProbs,structure,canPair, basepairWeighting);
			processOutsideSectors(sectors);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public ArrayList<Sector> createInsideSectors(int seqLen, InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs,
			int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting)
			{
		int divisions = (int)(threads*2);
		ArrayList<Sector> sectors = new ArrayList<Sector>();
		ParallelInfo info = new ParallelInfo(sectors);
		for (int b = 1; b <= seqLen; b += 1) {
			int j_divisions = Math.max(2, (seqLen-b)/divisions);			
			for (int j = 0; j <= seqLen-b; j += j_divisions) {
				Sector s = new Sector(sectors.size(), b, Math.min(seqLen, b+1), j,Math.min(seqLen-b, j+j_divisions), 0, seqLen, iProbs, pairingProbs,structure,canPair, basepairWeighting, info);
				sectors.add(s);
			}
		}

		return sectors;
			}

	public ArrayList<Sector> createOutsideSectors(int seqLen, InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
			NucleotideProbsPrecise pairingProbs, int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting)
			{
		int divisions = (int)(threads*2);
		ArrayList<Sector> sectors = new ArrayList<Sector>();


		ParallelInfo info = new ParallelInfo(sectors);
		for (int b = structure.length-2; b >= 0; b--) 
		{
			int j_divisions = Math.max(2, (seqLen-b)/divisions);
			for (int j = 0; j <= seqLen-b; j += j_divisions) {
				Sector s = new Sector(sectors.size(), b, Math.max(0, b-1), j,Math.min(seqLen-b, j+j_divisions), 0, seqLen, iProbs, oProbs, pairingProbs,structure,canPair, basepairWeighting, info);
				sectors.add(s);
			}
		}

		return sectors;
			}

	public void runInsideJob(Sector sector) {

		ParallelInfo info = sector.info;
		while(!info.canRunInside(sector))
		{

		}
		
		PointResUpperMatrix basepairWeighting = sector.basepairWeighting;
		if(basepairWeighting == null)
		{	
			InsideOutsideProbabilities iProbs = sector.iProbs;
			NucleotideProbsPrecise pairingProbs = sector.pairingProbs;
			boolean [][] canPair = sector.canPair;
			for (int b = sector.startb; b < sector.endb ; b++) {
				for (int j = sector.startj ; j < sector.endj ; j++) {
					//System.out.println(b+"\t"+j+"\t");
					// deal with contribution of rules of type U->VW
					for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
						PointRes tmp = PointRes.ZERO;
						for (int h = j ; h < j+b  ; h++) {
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
							probIncrement = rule2.getProbability()
									.multiply(pairingProbs.getPairingProbability(j, j+b))
									.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
							iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
						}
					}
				}	
			}
		}
		else
		{
			InsideOutsideProbabilities iProbs = sector.iProbs;
			NucleotideProbsPrecise pairingProbs = sector.pairingProbs;
			boolean [][] canPair = sector.canPair;
			for (int b = sector.startb; b < sector.endb ; b++) {
				for (int j = sector.startj ; j < sector.endj ; j++) {
					//System.out.println(b+"\t"+j+"\t");
					// deal with contribution of rules of type U->VW
					for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
						PointRes tmp = PointRes.ZERO;
						for (int h = j ; h < j+b  ; h++) {
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
							probIncrement = rule2.getProbability().multiply(basepairWeighting.get(j, j+b))
									.multiply(pairingProbs.getPairingProbability(j, j+b))
									.multiply(iProbs.getProb(rule2.getRight()[1],j+1,j+b-1));
							iProbs.increment(rule2.getLeft(), j, j+b, probIncrement);
						}
					}
				}	
			}
		}
	}	


	public void runOutsideJob(Sector sector) {

		
		ParallelInfo info = sector.info;
		while(!info.canRunOutside(sector))
		{

		}

		int seqLen = sector.structure.length;

		InsideOutsideProbabilities insideProbs = sector.iProbs;
		InsideOutsideProbabilities oProbs = sector.oProbs;
		NucleotideProbsPrecise pairingProbs = sector.pairingProbs;
		boolean [][] canPair = sector.canPair;

		int b = sector.startb;	
		PointResUpperMatrix basepairWeighting = sector.basepairWeighting;
		if(basepairWeighting == null)
		{		
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
						//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability()
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}
		else
		{
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
						//PointRes exp = PointRes.valueOf(Math.exp(-distances[j-1][j+b+1]/weight));
						probIncrement = rule2.getProbability().multiply(basepairWeighting.get(j-1, j+b+1))
								.multiply(pairingProbs.getPairingProbability(j-1, j+b+1))
								.multiply(oProbs.getProb(rule2.getLeft(), j-1, j+b+1));
						oProbs.increment(rule2.getRight()[1], j, j+b, probIncrement);
					}
				}
			}
		}

	}
	public List<Object> processInsideSectors(final List<Sector> inputs)
			throws InterruptedException, ExecutionException {

		ExecutorService service = Executors.newFixedThreadPool(threads);

		List<Future<Object>> futures = new ArrayList<Future<Object>>();
		for (final Sector input : inputs) {
			Callable<Object> callable = new Callable<Object>() {
				public Object call() throws Exception {
					input.running  = true;
					runInsideJob(input);
					input.finished = true;
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

	public List<Object> processOutsideSectors(final List<Sector> inputs)
			throws InterruptedException, ExecutionException {

		ExecutorService service = Executors.newFixedThreadPool(threads);

		List<Future<Object>> futures = new ArrayList<Future<Object>>();
		for (final Sector input : inputs) {
			Callable<Object> callable = new Callable<Object>() {
				public Object call() throws Exception {
					input.running  = true;
					runOutsideJob(input); 
					input.finished = true;
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


		public synchronized boolean canRunInside(Sector s)
		{
			Sector upto = getMostCompleteSector();
			if(s.startb == upto.startb || (upto.startb+1 >= s.startb && s.endj  < upto.startj) || s.startb-1 == 0 || upto.n + 1 == s.n)
			//if(s.startb == upto.startb || (upto.startb+1 >= s.startb && s.endj  < upto.endj) || s.startb-1 == 0)
			//if(s.startb == upto.startb || s.startb-1 == 0 || upto.n + 1 == s.n)
			{
				return true;
			}
			return false;
		}
		
		public synchronized boolean canRunOutside(Sector s)
		{
			Sector upto = getMostCompleteSector();
			//if(s.startb == upto.startb || (upto.startb+1 >= s.startb && s.endj  < upto.endj) || s.startb-1 == 0 || upto.n + 1 == s.n)
			//if(s.startb == upto.startb || (upto.startb+1 >= s.startb && s.endj  < upto.endj) || s.startb-1 == 0)
			if(s.startb == upto.startb  || (upto.startb <= s.startb+1 && s.endj  < upto.startj) || upto.n + 1 == s.n || (s.startb == 0 && s.startj == s.structure.length-1))
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
		final int [] structure;
		final boolean [][] canPair;
		final ParallelInfo info;
		final PointResUpperMatrix basepairWeighting;


		public Sector(int n, int startb, int endb, int startj, int endj, int starth, int endh, InsideOutsideProbabilities iProbs,
				NucleotideProbsPrecise pairingProbs,
				int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting, ParallelInfo info) {
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
			this.structure = structure;
			this.canPair = canPair;
			this.oProbs = null;
			this.info = info;
			this.basepairWeighting = basepairWeighting;
		}

		public Sector(int n, int startb, int endb, int startj, int endj, int starth, int endh, InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
				NucleotideProbsPrecise pairingProbs,
				int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting, ParallelInfo info) {
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
			this.structure = structure;
			this.canPair = canPair;
			this.info = info;
			this.basepairWeighting = basepairWeighting;
		}

		@Override
		public String toString() {
			return "Sector [startb=" + startb + ", endb=" + endb + ", startj="
					+ startj + ", endj=" + endj+"]";
		}
	}
	

/*
	public InsideOutsideProbabilities insideParallel (NucleotideProbsPrecise pairingProbs, double alpha, double tau, int [] structure, boolean [][] canPair)
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
	}*/

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.inoutside.IOsideCalculator#outside(uk.ac.ox.osscb.InsideOutsideProbabilities, uk.ac.ox.osscb.domain.NucleotideProbsPrecise, int[][], double, boolean[][])
	 */
	/*
	public InsideOutsideProbabilities outsideParallel(InsideOutsideProbabilities insideProbs, 
			NucleotideProbsPrecise pairingProbs, double alpha, double tau, int[] structure, boolean[][] canPair) {

		final int seqLen = structure.length;
		InsideOutsideProbabilities oProbs = new InsideOutsideProbabilities(grammar.getNonterminals(), seqLen);
		oProbs.setProb(this.grammar.getNonterminals()[0], 0, seqLen - 1, PointRes.ONE);
		runOutsideJobs(insideProbs, oProbs,pairingProbs,alpha,tau,structure,canPair);

		return oProbs;
	}*/



	/*
	public void runInsideJobs(InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{

		try {
			ArrayList<Sector> sectors = createInsideSectors(structure.length, iProbs, pairingProbs,alpha,tau,structure,canPair);
			processInsideSectors(sectors);

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	

	public void runInsideJobs(InsideOutsideProbabilities iProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair, PointResUpperMatrix basepairWeighting)
	{

		try {
			ArrayList<Sector> sectors = createInsideSectors(structure.length, iProbs, pairingProbs,alpha,tau,structure,canPair);
			processInsideSectors(sectors);

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	public void runOutsideJobs(InsideOutsideProbabilities iProbs, InsideOutsideProbabilities oProbs,
			NucleotideProbsPrecise pairingProbs, double alpha, double tau,
			int[] structure, boolean[][] canPair)
	{

		try {
			ArrayList<Sector> sectors = createOutsideSectors(structure.length, iProbs, oProbs, pairingProbs,alpha,tau,structure,canPair);
			processOutsideSectors(sectors);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}*/
}
