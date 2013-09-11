package uk.ac.ox.osscb;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ParallelInsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.ParallelValidatingIOSideCalculator;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

/**
 * Functional entry point. Does not parse input arguments, etc.
 * 
 * @author Vladimir, lepuslapis
 *
 */
public class KineticFold2 {
	
	private final Logger log = LoggerFactory.getLogger(KineticFold2.class);
	
	private boolean saveHelix = false; 
	
	public double sigma = 0; 
	
	public static void main(String [] args)
	{
		System.out.println(1/Double.POSITIVE_INFINITY);
	}


	/**
	 * 
	 * @param alignmentFile
	 * @param grammarFile
	 * @param paramsFile
	 * @param weight
	 */
	public void fold(String alignmentFile, String grammarFile, String paramsFile, double weight){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}
		
		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParser2().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile, parameters.getSAlphabet());
				
		// int[][] alignment = new AlignmentConverter().convert(align, parameters);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters); 

		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		
		Structure struct = new Structure(structure);
		
		//IterationsGenerator iterationsGenerator = new IterationsGenerator ();
		
		IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities completePPProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilities originalProbs = completePPProbs;
		double[][] distances = new DistancesCalculator2().distCalc(structure, completePPProbs);
		PosteriorProbabilities currentPostProbs = completePPProbs;
		
		
		
		boolean exitBecauseOfDiff = false;
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			canPair = new PossiblePairFinder().canPair(structure);				
			//completePPProbs = ppCalc.calculate(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
			completePPProbs = ppCalc.calculateParallel(insideProbs, outsideProbs, alignmentProbs, distances, weight, structure, canPair);
			
			PPOutputInternalResult2 postProbs = kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();
			
			//if (ppProbs.getDiff().signum()>0) {
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				//structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				structure = struct.makeNewStructure(ppProbs, true);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			OutputGenerator outputGenerator = new LoggingOutputGenerator();
			//outputGenerator.generate(structure);
			//dumpStructure(structure);
			outputGenerator.generate(struct);
			dumpStructure(struct);
		}
		

		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		//outputGenerator.generate(structure);
		//outputGenerator.generateFinal(structure);
		outputGenerator.generate(struct);
		outputGenerator.generateFinal(struct);
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			
			probs.savePosteriorProbabilities(new File(alignmentFile+".noevol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		writeDotBracketFile(new File(alignmentFile+".noevol.dbn"),new File(alignmentFile).getName(), struct);
	}
	
	public static void writeDotBracketFile(File dbnFile, String title, int [] structure)
	{
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(dbnFile));
			writer.write(">"+title+"\n");
			writer.write("\n");
			writer.write(LoggingOutputGenerator.dumpStructure(structure)+"\n");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void writeDotBracketFile(File dbnFile, String title, Structure struct)
	{
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(dbnFile));
			writer.write(">"+title+"\n");
			writer.write("\n");
			writer.write(LoggingOutputGenerator.dumpStructure(struct)+"\n");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	Random random= new Random(983118011408130222L);
	public String [] selectN(String [] align, int n)
	{
		int select = Math.min(align.length, n);
		boolean [] used = new boolean[align.length];
		String [] ret = new String[select];
		System.out.println("ret length "+ret.length+"\t"+n+"\t"+align.length);
		int count = 0;
		while(count < select)
		{
			int i = random.nextInt(align.length);
			if(!used[i])
			{
				System.out.println("getting sequence "+i);
				ret[count] = align[i];
				used[i] = true;
				count++;
			}
		}

		
		//align = ret;
		//boolean [] delete = new boolean[align[0].length()];
		/*Arrays.fill(delete, true);
		for(int i = 0 ; i < align.length ; i++)
		{
			String seq = align[i];
			for(int j = 0 ; j < seq.length() ; j++)
			{
				if(seq.charAt(j) != 'N' && seq.charAt(j) != '-')
				{
					delete[j] = false;
				}
			}
		}*/
		return ret;
	}
	
	public static double calculateReliablity(String [] inAlign, String alignmentFile, String grammarFile, String paramsFile, String treeFile)
	{
		boolean [] delete = CoFoldAnalogue.getGappyColumns(inAlign, 0.75);
		String [] align = CoFoldAnalogue.deleteColumns(inAlign, delete);
		
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		
		Grammar grammar = new GrammarParser().parse(grammarFile);
		IOsideCalculator ioCalc = new ParallelValidatingIOSideCalculator(new ParallelInsideOutsideCalculator(grammar), grammar);
		
		

		boolean[][] canPair = new PossiblePairFinder().canPair(structure, delete);
		
		
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		NucleotideProbsPrecise alignmentProbsEvol = EvolProbabilitiesCalculator.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile));
	
		InsideOutsideProbabilities insideProbsPPfoldAnalogue = ioCalc.insideE(alignmentProbsEvol, structure, canPair);
		InsideOutsideProbabilities outsideProbsPPfoldAnalogue = ioCalc.outsideE(insideProbsPPfoldAnalogue, alignmentProbsEvol, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalcPPfoldAnalogue = new PosteriorProbabilitiesCalculator(grammar);
		//int [] structuree
		int[] structurePPfoldAnalogue = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structurePPfoldAnalogue.length; posIdx++){
			structurePPfoldAnalogue[posIdx] = Constants.UnpairedBaseIdx;
		}
		PosteriorProbabilities currentPostProbsAnalogue = ppCalcPPfoldAnalogue.calculateE(insideProbsPPfoldAnalogue, outsideProbsPPfoldAnalogue, alignmentProbsEvol, structurePPfoldAnalogue, canPair);
		
		
		double [][] basePairProbabilities = currentPostProbsAnalogue.getBasePairProbs();
		double [] single = StructureData.getSingleProbs(basePairProbabilities);
		double reliability = 0;
		structurePPfoldAnalogue = RNAFoldingTools.getPosteriorDecodingConsensusStructure(basePairProbabilities, single);
		for(int k = 0 ; k < structurePPfoldAnalogue.length ; k++)
		{
			if(structurePPfoldAnalogue[k] == 0)
			{
				reliability += single[k];
			}
			else
			{
				reliability += basePairProbabilities[k][structurePPfoldAnalogue[k]-1];
			}
		}
		
		reliability /= (double)(structurePPfoldAnalogue.length);
		
		return reliability;
	}
	
	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight, double delta2, double p){
		weight = Constants.weight; 
		double gap_perc = 0.75;				
		
		int numWeak = 0; 
		int numStrong = 0; 
		int numOx = 0; 
		
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}
		


		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);

		double ppfoldReliablity = calculateReliablity(align, alignmentFile, grammarFile, paramsFile, treeFile);
		

		//align = selectN(align, (int)delta2);

		//boolean [] delete =new boolean[align[0].length()];
		boolean [] delete = CoFoldAnalogue.getGappyColumns(align, gap_perc);
		align = CoFoldAnalogue.deleteColumns(align, delete);
		//System.out.println(align[0]);
		//
		//.out.println(align[0]);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		//NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);

		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		NucleotideProbsPrecise alignmentProbsEvol = EvolProbabilitiesCalculator2.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile), p);
		
		

		double sigma = Constants.sigma; 
		sigma = 0.01; 
		boolean cotranscriptional = true;
		if(ppfoldReliablity < 0.7)
		{
			/*weight = Double.POSITIVE_INFINITY;
			Constants.IterationCutOffDf = PointRes.valueOf(-0.2);
			Constants.IterationCutOffDouble = 0.0;
			Constants.IterationCutOffDr = PointRes.valueOf(0.15);
			*/
			weight = Double.POSITIVE_INFINITY;
			Constants.IterationCutOffDf = PointRes.valueOf(-0.2);
			Constants.IterationCutOffDouble = 0.0;
			Constants.IterationCutOffDr = PointRes.valueOf(0.15);
			System.out.println("hit 0.7");
		}
		else
			if(ppfoldReliablity < 0.85)
			{
				/*weight = 0.5;
				Constants.IterationCutOffDf = PointRes.valueOf(0.1);
				Constants.IterationCutOffDouble = 0.4;
				Constants.IterationCutOffDr = PointRes.valueOf(0.7);
				*/
				weight = 0.5;
				Constants.IterationCutOffDf = PointRes.valueOf(0.4);
				Constants.IterationCutOffDouble = 0.6;
				Constants.IterationCutOffDr = PointRes.valueOf(0.9);
				System.out.println("hit 0.85");
			}
			else
		if(ppfoldReliablity < 0.90)
		{
			/*weight = 0.5;
			Constants.IterationCutOffDf = PointRes.valueOf(0.1);
			Constants.IterationCutOffDouble = 0.4;
			Constants.IterationCutOffDr = PointRes.valueOf(0.7);
			*/
			weight = 0.5;
			//sigma = 0.01; 
			Constants.IterationCutOffDf = PointRes.valueOf(0.5);
			Constants.IterationCutOffDouble = 0.6;
			Constants.IterationCutOffDr = PointRes.valueOf(1.0);
			System.out.println("hit 0.9");
		}
		else
		if(ppfoldReliablity < 0.95)
		{
			/*weight = 0.5;
			Constants.IterationCutOffDf = PointRes.valueOf(0.1);
			Constants.IterationCutOffDouble = 0.4;
			Constants.IterationCutOffDr = PointRes.valueOf(0.7);
			*/
			weight = Constants.weight;
			Constants.IterationCutOffDf = PointRes.valueOf(0.5);
			Constants.IterationCutOffDouble = 0.6;
			Constants.IterationCutOffDr = PointRes.valueOf(1.0);
			System.out.println("hit 0.95");
		}
		else
		if(ppfoldReliablity < 0.99)
		{
			weight = Constants.weight;
			Constants.IterationCutOffDf = PointRes.valueOf(0.5);
			Constants.IterationCutOffDouble = 0.6;
			Constants.IterationCutOffDr = PointRes.valueOf(1.0);
			System.out.println("hit 0.99");
		}
		else
		{
			weight = Constants.weight;
			Constants.IterationCutOffDf = PointRes.valueOf(0.5);
			Constants.IterationCutOffDouble = 0.6;
			Constants.IterationCutOffDr = PointRes.valueOf(1.0);
			// should do ppfold here
			System.out.println("hit else");
		}
		Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble); // important to do this assignment
		//Constants.IterationCutOffDr = Constants.IterationCutOffDf;
		
		if(log.isDebugEnabled()){
			log.debug(String.format("Pairing probs (NON-evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (NON-evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
			
			log.debug(String.format("Pairing probs (Evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (Evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
		}
		
		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		int[] structure2 = structure; 
		Structure struct = new Structure(structure);
		
		//IterationsGenerator iterationsGenerator = new IterationsGenerator ();
		//
		//IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		IOsideCalculator ioCalc = new ParallelValidatingIOSideCalculator(new ParallelInsideOutsideCalculator(grammar), grammar);
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);
		
		
		boolean evol = false;
		//boolean evol = true; 
		boolean exitBecauseOfDiff = false;
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure, delete);
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		int fastIterations = 0;
		int iterSoFar = 0;

		int jstar = 0;
		ArrayList<Pair> stopPair = new ArrayList<Pair>();
		for(; iterSoFar < Constants.MaxIterations; iterSoFar++){
			//canPair = new PossiblePairFinder().canPair(structure);
			canPair = new PossiblePairFinder().canPair(structure, delete);
			for(Pair pair : stopPair)
			{
				if(pair.i == -1)
				{
					for(int i = 0 ; i < canPair.length ; i++)
					{
						canPair[i][pair.j] = false;
						canPair[pair.j][i] = false;
					}
				}
			}
					
			PPOutputInternalResult2 postProbs = evol ? 
					kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();

			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			
			
			PointRes[] unpairedProbs = currentPostProbs.getUnpairedProbs();
			PointRes[][] pairedProbs = currentPostProbs.getPairedProbs();
			PointRes[][] diffs = ppCalc.getDiffs(pairedProbs, unpairedProbs, canPair);

			Constants.mu2 = PointRes.valueOf(0);

			Constants.currentIterationCutOff = Constants.IterationCutOff;

			PPOutput oxfoldPProbs = ppProbs;
			
			boolean madeBase = true;
			this.saveHelix = false; 
			if(cotranscriptional)
			{
				System.out.println("Co-transcriptional foldevol");
				//set dr df
				//ArrayList<BasePair> pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOff);
				ArrayList<BasePair> pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOffDr);
				//this.saveHelix = true; 
				System.out.println(pairsList);
				BasePair next = findNextCotranscriptionalBasePair(pairsList,jstar, sigma);
				//Constants.currentIterationCutOff = Constants.IterationCutOffDr; 
				 		
				//BasePair next = findNextCotranscriptionalBasePair(canPair, pairedProbs, diffs, Constants.IterationCutOff, sigma, jstar);
				if(next != null)
				{
					madeBase = true; 
					this.saveHelix = true;
					Constants.currentIterationCutOff = Constants.IterationCutOffDr; 
					boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, next.i, next.j);
					Helix helix = new HelicesMaker().makeHelix(next.i,next.j,diffs,canPair);
					PointRes diff = diffs[next.i][next.j];
					PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,next.i,next.j,pairedProbs);
				
					System.out.println("Next "+jstar+"\t"+next);
					ppProbs = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
					numStrong++; 
					jstar = Math.max(jstar, helix.getRightIdx());
				}
				else{
					//look for df satisfied bases
					System.out.println("looking for weak helices");
					//Constants.currentIterationCutOff = Constants.IterationCutOffDf; 
					pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOffDf);
					this.saveHelix = false; 
					System.out.println(pairsList);
					next = findNextCotranscriptionalBasePair(pairsList,jstar, sigma);
							
					//BasePair next = findNextCotranscriptionalBasePair(canPair, pairedProbs, diffs, Constants.IterationCutOff, sigma, jstar);
					if(next != null)
					{
						madeBase = true; 
						Constants.currentIterationCutOff = Constants.IterationCutOffDf;
						Constants.mu2 = PointRes.valueOf(Math.min(0.0, Constants.IterationCutOffDf.toDouble()));
						System.out.println("looking for weak helices: dealing with next basepair");
						boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, next.i, next.j);
						Helix helix = new HelicesMaker().makeHelix(next.i,next.j,diffs,canPair);
						PointRes diff = diffs[next.i][next.j];
						PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,next.i,next.j,pairedProbs);
					
						System.out.println("Next "+jstar+"\t"+next);
						ppProbs = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
						numWeak++; 
						jstar = Math.max(jstar, helix.getRightIdx());
					}
				}
			}
			
			if (madeBase && ppProbs.getDiff().compareTo(Constants.currentIterationCutOff)>0 && oxfoldPProbs.getDiff().compareTo(Constants.currentIterationCutOff)>0) {				
			//if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0 && oxfoldPProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {				
				//structure = new StructureUtils().makeNewStructure(structure, ppProbs);				
				structure = struct.makeNewStructure(ppProbs, this.saveHelix);		
				System.out.println("added to structure");
				System.out.println("current delta: " + Constants.currentIterationCutOff);
				System.out.println("ppProbs.getDiff(): "+ ppProbs.getDiff());
				double currentDelta = ppProbs.getDiff().doubleValue();
				
			} else {
				//global step
				this.saveHelix = true; 
				//get rid of weak helices
				Constants.currentIterationCutOff = Constants.IterationCutOff; 
				
				int[] oldStructure = struct.getPairings(); 
				System.out.print("old structure: ");
				for (int a = 0; a < oldStructure.length; a++){
					System.out.print(oldStructure[a] + "\t");
				}
				System.out.println();
				structure = struct.removeWeakPairs(); 
				
				System.out.print("new structure: ");
				for (int a = 0; a < structure.length; a++){
					System.out.print(structure[a] + "\t");
				}
				System.out.println(); 
				
				
				canPair = new PossiblePairFinder().canPair(structure, delete);


				InsideOutsideProbabilities oldProbs = insideProbs;
				try{
					insideProbs = ioCalc.insideE(alignmentProbs, oldStructure, canPair);
				}
				catch(Exception e){
					insideProbs = oldProbs; 
					System.err.println("inside prob calc FAILS: " + structure.toString());
					
					try {
						BufferedWriter writerFinal = new BufferedWriter(new FileWriter( "outputAOxfold/insideFails.txt", true));
						writerFinal.write(alignmentFile + "\t" +LoggingOutputGenerator.dumpStructure(structure)+"\n");
						writerFinal.close();
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					
				}
				
				outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, oldStructure, canPair);
				ppCalc = new PosteriorProbabilitiesCalculator(grammar);
				currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, oldStructure, canPair);
				boolean first = true; 
				
				for(; iterSoFar < Constants.MaxIterations; iterSoFar++){
					//canPair = new PossiblePairFinder().canPair(structure);
					canPair = new PossiblePairFinder().canPair(structure, delete);
					for(Pair pair : stopPair)
					{
						if(pair.i == -1)
						{
							for(int i = 0 ; i < canPair.length ; i++)
							{
								canPair[i][pair.j] = false;
								canPair[pair.j][i] = false;
							}
						}
					}
						
					if(first){
						postProbs = evol ? 
								kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, oldStructure, currentPostProbs) 
									:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, oldStructure, currentPostProbs);
						
						ppProbs = postProbs.getPpProbs();
						currentPostProbs = postProbs.getPosteriorProbs();
						first = false; 
					}
					else{
						postProbs = evol ? 
								kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
									:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
						
						ppProbs = postProbs.getPpProbs();
						currentPostProbs = postProbs.getPosteriorProbs();
					}
					

					if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
						evol = true;
						continue;
					}
					
					
					unpairedProbs = currentPostProbs.getUnpairedProbs();
					pairedProbs = currentPostProbs.getPairedProbs();
					diffs = ppCalc.getDiffs(pairedProbs, unpairedProbs, canPair);

					oxfoldPProbs = ppProbs;
					
					System.out.println("ppProbs diff: "+ ppProbs.getDiff().toDouble());
					System.out.println("oxfoldPProbs diff: "+ oxfoldPProbs.getDiff().toDouble());
					if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0 && oxfoldPProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {				
						//structure = new StructureUtils().makeNewStructure(structure, ppProbs);		
						System.out.println("Global: new structure made");
						structure = struct.makeNewStructure(ppProbs, this.saveHelix);
						numOx++; 
						dumpCurrentOutput(ppProbs);
						System.out.println("added to structure g");
					} else {
						//if(structure.equals(structure2)){
							exitBecauseOfDiff = true;
							dumpCurrentOutput(ppProbs);
							break;
						/*}
						else{
							structure = structure2; 
						}*/
					}
				}
				
				//exitBecauseOfDiff = true;
				//dumpCurrentOutput(ppProbs);
				//if(structure.equals(structure2)){
					break;
				//}
			}
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			dumpStructure(struct);
		}



		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		

		structure = CoFoldAnalogue.reinsertDeleted(structure, delete, Constants.UnpairedBaseIdx);
		struct.setPairings(structure);
		outputGenerator.generate(struct);
		outputGenerator.generateFinal(struct);
		
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			// Added in
			probs.pairedProbs = CoFoldAnalogue.reinsertDeleted(currentPostProbs.pairedProbs, delete, PointRes.ZERO);
			structure = CoFoldAnalogue.reinsertDeleted(structure, delete, Constants.UnpairedBaseIdx);
			probs.savePosteriorProbabilities(new File(alignmentFile+".evol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String metadata = ";oxfolditer="+iterSoFar+",fastiter="+fastIterations+",delta2="+delta2;
		System.out.println(metadata);
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName()+metadata, struct);
		System.out.println("Strong helix = " + numStrong);
		System.out.println("Weak helix = " + numWeak);
		System.out.println("OxfoldI helix = " + numOx);
	}

	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight, double delta2){
		//weight = Double.POSITIVE_INFINITY;
		//weight = 0.5;
		weight = Constants.weight; 
		//	double sigma = delta2;
		double gap_perc = 2;
		//double sigma = 0.1;
		double sigma = Constants.sigma; 
		boolean cotranscriptional = false;
		
		boolean multipleHelixFormation = false;
		
		Constants.mu2 = PointRes.valueOf(0);
		
		
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}
		


		Grammar grammar = new GrammarParser().parse(grammarFile);
		
		EvolutionaryParameters parameters = new ParameterParserEvolutionary().parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		

		//align = selectN(align, (int)delta2);

		//boolean [] delete =new boolean[align[0].length()];
		boolean [] delete = CoFoldAnalogue.getGappyColumns(align, gap_perc);
		align = CoFoldAnalogue.deleteColumns(align, delete);
		//System.out.println(align[0]);
		//
		//.out.println(align[0]);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		//NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);

		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		NucleotideProbsPrecise alignmentProbsEvol = EvolProbabilitiesCalculator.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile));
		
		
		
		// Added in
		//
		
		
		//ArrayList<String> sequences = new ArrayList<String>();
		//ArrayList<String> sequenceNames = new ArrayList<String>();
		//IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		//alignmentProbsEvol = EvolProbabilitiesCalculator.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile));
		//alignmentProbsEvol = alignmentProbs;
		
		
		if(log.isDebugEnabled()){
			log.debug(String.format("Pairing probs (NON-evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (NON-evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
			
			log.debug(String.format("Pairing probs (Evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (Evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
		}
		
		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		
		//Structure struct = new Structure(structure);
		
		//IterationsGenerator iterationsGenerator = new IterationsGenerator ();
		//
		//IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		IOsideCalculator ioCalc = new ParallelValidatingIOSideCalculator(new ParallelInsideOutsideCalculator(grammar), grammar);
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);
		
		
		boolean evol = false;
		boolean exitBecauseOfDiff = false;
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure, delete);
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		int fastIterations = 0;
		int iterSoFar = 0;

		int jstar = 0;
		ArrayList<Pair> stopPair = new ArrayList<Pair>();
		for(; iterSoFar < Constants.MaxIterations; iterSoFar++){
			//canPair = new PossiblePairFinder().canPair(structure);
			canPair = new PossiblePairFinder().canPair(structure, delete);
			for(Pair pair : stopPair)
			{
				if(pair.i == -1)
				{
					for(int i = 0 ; i < canPair.length ; i++)
					{
						canPair[i][pair.j] = false;
						canPair[pair.j][i] = false;
					}
				}
			}
					
			PPOutputInternalResult2 postProbs = evol ? 
					kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();

			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			
			
			PointRes[] unpairedProbs = currentPostProbs.getUnpairedProbs();
			PointRes[][] pairedProbs = currentPostProbs.getPairedProbs();
			PointRes[][] diffs = ppCalc.getDiffs(pairedProbs, unpairedProbs, canPair);

			PPOutput oxfoldPProbs = ppProbs;
			if(cotranscriptional)
			{
				System.out.println("Co-transcriptional");
				ArrayList<BasePair> pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOff);
				
				BasePair next = findNextCotranscriptionalBasePair(pairsList,jstar, sigma);
						
				//BasePair next = findNextCotranscriptionalBasePair(canPair, pairedProbs, diffs, Constants.IterationCutOff, sigma, jstar);
				if(next != null)
				{
					boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, next.i, next.j);
					Helix helix = new HelicesMaker().makeHelix(next.i,next.j,diffs,canPair);
					PointRes diff = diffs[next.i][next.j];
					PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,next.i,next.j,pairedProbs);
				
					System.out.println("Next "+jstar+"\t"+next);
					ppProbs = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
					
					jstar = Math.max(jstar, helix.getRightIdx());
				}
			}
			
			
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0 && oxfoldPProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {				
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);				
				//structure = struct.makeNewStructure(ppProbs, false);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			//dumpStructure(struct);
			dumpStructure(structure);
		}
		

		
		// remove low confidence helices, run again

		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		//outputGenerator.generate(struct);
		//outputGenerator.generateFinal(struct);
		
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			// Added in
			probs.pairedProbs = CoFoldAnalogue.reinsertDeleted(currentPostProbs.pairedProbs, delete, PointRes.ZERO);
			structure = CoFoldAnalogue.reinsertDeleted(structure, delete, Constants.UnpairedBaseIdx);
			probs.savePosteriorProbabilities(new File(alignmentFile+".evol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String metadata = ";oxfolditer="+iterSoFar+",fastiter="+fastIterations+",delta2="+delta2;
		System.out.println(metadata);
		//writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName()+metadata, struct);
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName()+metadata, structure);
	}

	public static ArrayList<BasePair> listPossibleBasePairs(boolean [][] canPair, PointRes [][] pairedProbs, PointRes [][] diffs, PointRes deltaGreaterThan)
	{
		ArrayList<BasePair> pairs = new ArrayList<BasePair>();
		
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = i + 1 ; j < canPair.length ; j++)
			{
				if(canPair[i][j] && diffs[i][j].compareTo(deltaGreaterThan) > 0)
				{
					pairs.add(new BasePair(i, j, pairedProbs[i][j], diffs[i][j]));
				}
			}
		}
		
		return pairs;
	}
	

	
	public static BasePair findNextCotranscriptionalBasePair(ArrayList<BasePair> possibleBasePairs, int jstar, double kfactor)
	{
		PointRes k = new PointRes(kfactor);
		
		BasePair minPair = null;
		PointRes min = null;

		System.out.println("SIGMA = "+k);
		for(BasePair basePair : possibleBasePairs)
		{
			PointRes j =new PointRes(Math.max(basePair.i, basePair.j));
			//PointRes cotranscriptional = k.multiply(j);
			PointRes cotranscriptional = k.multiply(PointRes.max(PointRes.ZERO, j.subtract(new PointRes(jstar)).subtract(Constants.AYOffset)));
			//add -25 here
			if(basePair.probability.compareTo(PointRes.ZERO) == 0)
			{
				continue;
			}
			PointRes val = PointRes.ONE.divide(basePair.probability).add(cotranscriptional);
			//PointRes val = PointRes.ONE.divide(basePair.probability).add(k.multiply(cotranscriptional));
			if(min == null || val.compareTo(min) < 0)
			{
				System.out.println("COT = "+PointRes.ONE.divide(basePair.probability)+"\t"+cotranscriptional);
				minPair = basePair;
				min = val;
			}
		}
		return minPair;
	}
	
	public static BasePair findNextCotranscriptionalBasePair(boolean [][] canPair, PointRes [][] pairedProbs, PointRes [][] diffs, PointRes deltaGreaterThan, double kfactor, int jstar)
	{
		PointRes k = new PointRes(kfactor);
		
		BasePair minPair = null;
		PointRes min = null;
		
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = i + 1 ; j < canPair.length ; j++)
			{
				if(canPair[i][j] && diffs[i][j].compareTo(deltaGreaterThan) > 0)
				{
					BasePair temp = new BasePair(i, j, pairedProbs[i][j], diffs[i][j]);
					PointRes jval = new PointRes(Math.max(temp.i, temp.j));
					//PointRes cotranscriptional = k.multiply(j);
					PointRes cotranscriptional = k.multiply(PointRes.max(PointRes.ZERO, jval.subtract(new PointRes(jstar))));
					
					if(temp.probability.compareTo(PointRes.ZERO) == 0)
					{
						continue;
					}
					PointRes val = PointRes.ONE.divide(temp.probability).add(cotranscriptional);
					//PointRes val = PointRes.ONE.divide(temp.probability).add(k.multiply(cotranscriptional));
					if(min == null || val.compareTo(min) < 0)
					{
						minPair = temp;
						min = val;
					}
				}
			}
		}

		return minPair;
	}
	
	public static class Pair
	{
		public int i;
		public int j;
		
		public Pair(int i, int j)
		{
			this.i = i;
			this.j = j;
		}
	}
	
	public static class BasePair
	{
		public int i;
		public int j;
		PointRes probability;
		PointRes delta;
		
		public BasePair(int i, int j, PointRes probability, PointRes delta)
		{
			this.i = i;
			this.j = j;
			this.probability = probability;
			this.delta = delta;
		}

		@Override
		public String toString() {
			return "BasePair [i=" + i + ", j=" + j + ", probability="
					+ probability + ", delta=" + delta + "]";
		}
		
		
	}
	
	public static PPOutput getMaxBasePair(boolean [][] canPair, PointRes [][] pairedProbs, int [] structure)
	{
		int leftIdx = -1;
		int rightIdx = -1;
		
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = i + 1 ; j < canPair.length ; j++)
			{
				if(canPair[i][j]  && structure[i] < 0 && structure[j] < 0)
				{
					if((leftIdx == -1 && rightIdx == -1) || pairedProbs[i][j].compareTo( pairedProbs[leftIdx][rightIdx]) > 0)
					{
						leftIdx = i;
						rightIdx = j;
					}
				}
			}
		}
		
		return new PPOutput(leftIdx, rightIdx, 0, new PointRes(0.0), new PointRes(0.0));
	}
	
	private void dumpExitReason(boolean exitBecauseOfDiff) {
		ProgramOutput.outMsg(String.format("\texiting because of %s", 
				exitBecauseOfDiff ? "difference under threshold." : "number of iterations."));
	}

	private void dumpCurrentOutput(PPOutput ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g\trprob: %g"
				, ppProbs.getLeftIdx()
				, ppProbs.getRightIdx()
				, ppProbs.gethelixLength()
				, ppProbs.getDiff().doubleValue()
				, ppProbs.getComp().doubleValue()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}
	
	private void dumpStructure(Structure structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}
}
