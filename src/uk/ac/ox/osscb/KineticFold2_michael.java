package uk.ac.ox.osscb;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.analysis.BasePairMetrics;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.RankingAnalyses;
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
public class KineticFold2_michael {
	
	private final Logger log = LoggerFactory.getLogger(KineticFold2_michael.class);
	
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
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			OutputGenerator outputGenerator = new LoggingOutputGenerator();
			outputGenerator.generate(structure);
			dumpStructure(structure);
		}
		

		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		try {
			PosteriorProbabilities probs = new PosteriorProbabilities(structure);
			
			probs.savePosteriorProbabilities(new File(alignmentFile+".noevol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		writeDotBracketFile(new File(alignmentFile+".noevol.dbn"),new File(alignmentFile).getName(), structure);
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
	
	public static double calculateDelta(double reliability)
	{
		double ret =0;
		if(reliability < 0.75)
		{
			ret = 6.4325*reliability - 4.2606;
			ret = Math.max(-0.4, ret);
			ret = Math.min(0, ret);
		}
		
		System.out.println("Reliability="+reliability+"\t"+",Delta="+ret);
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

	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double weight, double delta2){
		
		//double weight2 = weight;
		/*weight = Double.POSITIVE_INFINITY;
		//weight = 0.5;
		Constants.IterationCutOffDouble = 0.0;
		Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble);

		//	double sigma = delta2;
		double gap_perc = 0.25;
		double sigma = 0.1;
		boolean cotranscriptional = false;
		boolean multipleHelixFormation = false;
		*/


		//weight = Double.POSITIVE_INFINITY;
		weight = 0.5;
		//weight = 0.5;
		Constants.mu = 0.0;
		//Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble);
		//double gap_perc = 0.75;
		//double sigma = 0.05;
		double sigma = 0;
		boolean cotranscriptional = false;
		boolean multipleHelixFormation = false;
		System.out.println("RUNNING "+Constants.IterationCutOffDouble+"\t"+weight+"\t"+Constants.mu+"\t"+Constants.gapPercentage);

		
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
		boolean [] delete = CoFoldAnalogue.getGappyColumns(align, Constants.gapPercentage);
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

		//IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new InsideOutsideCalculator(grammar), grammar);
		IOsideCalculator ioCalc = new ParallelValidatingIOSideCalculator(new ParallelInsideOutsideCalculator(grammar), grammar);
	
		
		
		boolean evol = false;
		boolean exitBecauseOfDiff = false;
		
		boolean[][] canPair = new PossiblePairFinder().canPair(structure, delete);
		//alignmentProbs = alignmentProbsEvol;
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		
		
		/*
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
		System.out.println("ALIGNMENT SIZE="+align.length);
		reliability /= (double)(structurePPfoldAnalogue.length);
		System.out.println(alignmentFile);
		//Constants.IterationCutOffDouble = -0.3;
		//Constants.IterationCutOffDouble = calculateDelta(reliability);
		//Constants.IterationCutOffDouble = calculateDelta(0.0);
		//Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble);
		
		
		
		double [][] meaDiffs = PosteriorProbabilitiesCalculator.getDiffs(basePairProbabilities,canPair);
		int [] meaStructure = RNAFoldingTools.getPosteriorDecodingConsensusStructure(basePairProbabilities);
		ArrayList<Double> diffList = new ArrayList<Double>();
		double mindiff = Double.MAX_VALUE;
		double maxdiff = Double.MIN_VALUE;
		double minmaxdelta= 0.95;

		ArrayList<Helix> helices = SimonsIdea.getHelices(meaStructure, meaDiffs);
		for(Helix helix : helices)
		{
			if(helix.getHelixLength() >= 4)
			{
				minmaxdelta = Math.min(minmaxdelta, helix.maxdelta);
			}
		}
		
		for(int i = 0 ; i < meaStructure.length ; i++)
		{
			if(meaStructure[i] != 0 && i < meaStructure[i]-1)
			{
				diffList.add(meaDiffs[i][meaStructure[i]-1]);
				mindiff = Math.min(meaDiffs[i][meaStructure[i]-1], mindiff);
				maxdiff = Math.max(meaDiffs[i][meaStructure[i]-1], maxdiff);
			}*/
			
			
			/*if(meaStructure[i] != 0 && (i == 0 || meaStructure[i-1] == 0) && i < meaStructure[i]-1) // find first base-pair in helix
			{
				int index = i;
				double maxPij = Double.MIN_VALUE;
				double bestDiffInHelix = 0;
				while(index <= (meaStructure.length+1)/2 && meaStructure[index] != 0) // find base-pair in helix with maxij
				{
					if(basePairProbabilities[index][meaStructure[index]-1] > maxPij)
					{
						maxPij = basePairProbabilities[index][meaStructure[i]-1];
						bestDiffInHelix = meaDiffs[index][meaStructure[index]-1];
					}
					index++;
				}
				diffList.add(bestDiffInHelix);
			}*/
		/*}
		
		
		System.out.println("Differences="+mindiff+"\t"+maxdiff+"\t"+minmaxdelta);
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(alignmentFile+".deltas"));
			writer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(meaStructure));
			writer.newLine();
			for(int i = 0 ; i < diffList.size() ; i++)
			{
				writer.write(diffList.get(i)+"");
				if(i != meaDiffs.length-1)
				{
					writer.write(",");
				}
			}
			writer.newLine();
			for(int i = 0 ; i < meaDiffs.length ; i++)
			{
				for(int j = 0 ; j < meaDiffs.length ; j++)
				{
					writer.write(meaDiffs[i][j]+"");
					if(j != meaDiffs.length-1)
					{
						writer.write(",");
					}
				}
				writer.newLine();
			}
			writer.close();
			
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
		
		System.out.println("Reliability = "+reliability+"\t"+ppfoldReliablity);
		*/
		System.out.println("Reliability = "+ppfoldReliablity);
		if(ppfoldReliablity < 0.7)
		{
			weight = Double.POSITIVE_INFINITY;
			Constants.IterationCutOffDouble = 0.0;
		}
		else
		if(ppfoldReliablity < 0.90)
		{
			weight = 0.5;
			Constants.IterationCutOffDouble = 0.4;
		}
		else
		if(ppfoldReliablity < 0.99)
		{
			weight = 0.5;
			Constants.IterationCutOffDouble = 0.8;
		}
		else
		{
			weight = 0.5;
			Constants.IterationCutOffDouble = 0.8;
			// should do ppfold here
		}
		
		//weight = weight2;
	//	Constants.IterationCutOffDouble = 0.5;
		//System.out.println("Weight = "+weight+"\t"+Constants.IterationCutOffDouble);
		
		//weight = 0.5;
		//Constants.IterationCutOffDouble = 0.5;
	    Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble);
		
		
		KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);
		
		int fastIterations = 0;
		int iterSoFar = 0;

		int jstar = 0;
		for(; iterSoFar < Constants.MaxIterations; iterSoFar++){
			canPair = new PossiblePairFinder().canPair(structure, delete);
			
			PPOutputInternalResult2 postProbs = evol ? 
					kFPppCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						:kFPppCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();
			
/*
			
			try {
				File experimentalFile = new File(alignmentFile.replaceAll("_oxfold\\.fas$", "")+".experimental");
				BufferedReader buffer = new BufferedReader(new FileReader(experimentalFile));
				String exp = buffer.readLine();
				int [] sites = RNAFoldingTools.getPairedSitesFromDotBracketString(exp);
				buffer.close();

				int [] pred = new int[structure.length];
				for(int i = 0 ; i < pred.length ; i++)
				{
					pred[i] = structure[i]+1;
				}

				File deltaFile = new File(alignmentFile+".deltas.txt"); 
				if(iterSoFar == 0)
				{
					System.out.println("NEW FILE");
					IO.writeLine(deltaFile, "Delta\tSensitivity\tPPV\tFScore", true, false);
				}
				System.out.println("WRITE LINE "+iterSoFar);
				IO.writeLine(deltaFile, ""+ppProbs.getDiff().doubleValue()+"\t"+BasePairMetrics.calculateSensitivity(sites, pred)+"\t"+BasePairMetrics.calculatePPV(sites, pred)+"\t"+BasePairMetrics.calculateFScore(sites, pred), true, true);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			*/

			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			
			
			PointRes[] unpairedProbs = currentPostProbs.getUnpairedProbs();
			PointRes[][] pairedProbs = currentPostProbs.getPairedProbs();
			PointRes[][] diffs = ppCalc.getDiffs(pairedProbs, unpairedProbs, canPair);
			

			/*try {
				if(next != null)
				{
					IO.writeLine(new File("basepairs.txt"), next.toString(), true, true);
					IO.writeLine(new File("basepairs.txt"), ppProbs.toString()+"\n", true, true);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}*/
			

			
			PPOutput oxfoldPProbs = ppProbs;
			//BasePair maxPij = null;
			if(cotranscriptional)
			{
				System.out.println("Co-transcriptional");
				ArrayList<BasePair> pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOff);
				System.out.println(pairsList);
				BasePair next = findNextCotranscriptionalBasePair(pairsList,jstar, sigma);
						
				//BasePair next = findNextCotranscriptionalBasePair(canPair, pairedProbs, diffs, Constants.IterationCutOff, sigma, jstar);
				if(next != null)
				{
					boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, next.i, next.j);
					//Helix helix = new HelicesMaker().makeHelix(next.i,next.j,diffs,canPair);
					Helix helix = new HelicesMaker().makeHelix(next.i,next.j, pairedProbs,diffs,canPair);
					PointRes diff = diffs[next.i][next.j];
					PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,next.i,next.j,pairedProbs);
				
					System.out.println("Next "+jstar+"\t"+next);
					ppProbs = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
					
					jstar = Math.max(jstar, helix.getRightIdx());
				}
			}
			/*else
			{
				maxPij = findNextMaxPijBasePair(canPair, pairedProbs, diffs, Constants.IterationCutOff);
				if(maxPij != null)
				{
					boolean[][] incomp = new IncompatiblePairsFinder().find(canPair, maxPij.i, maxPij.j);
					Helix helix = new HelicesMaker().makeHelix(maxPij.i,maxPij.j,diffs,canPair);
					PointRes diff = diffs[maxPij.i][maxPij.j];
					PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,maxPij.i,maxPij.j,pairedProbs);
				
					ppProbs = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);

				}			
			}*/
		
			
			//System.out.println("Cutoff = "+Constants.IterationCutOffDouble+"\t"+ppProbs.getDiff().doubleValue());
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0 && oxfoldPProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				//Constants.IterationCutOffDouble= ppProbs.getDiff().doubleValue()*0.5;
				//Constants.IterationCutOffDouble = Math.max(Constants.IterationCutOffDouble, Math.exp(ppProbs.getDiff().doubleValue()-0.7) - 1);
				//Constants.IterationCutOffDouble = Math.max(Constants.IterationCutOffDouble, ppProbs.getDiff().doubleValue()*0.5);
			///	Constants.IterationCutOff = PointRes.valueOf(Constants.IterationCutOffDouble);		
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);	
				
				double currentDelta = ppProbs.getDiff().doubleValue();
				if(multipleHelixFormation)
				{
					for(int c = 0 ;  ; c++)
					{
						fastIterations++;
						if(ppProbs.getDiff().compareTo(new PointRes(delta2))<0)
						{
							break;
						}
						boolean [][] canPair2 = new PossiblePairFinder().canPair(structure);	
						/*for(Pair pair : stopPair)
						{
							if(pair.i == -1)
							{
								for(int i = 0 ; i < canPair2.length ; i++)
								{
									canPair2[i][pair.j] = false;
									canPair2[pair.j][i] = false;
								}
							}
						}*/
						unpairedProbs = currentPostProbs.getUnpairedProbs();
						pairedProbs = currentPostProbs.getPairedProbs();
						
	
						PPOutput nextBasePair = null;
						
						if(cotranscriptional)
						{
							/*ArrayList<BasePair> pairsList = listPossibleBasePairs(canPair, pairedProbs, diffs, Constants.IterationCutOff);
							System.out.println(pairsList);
							BasePair next = findNextCotranscriptionalBasePair(pairsList,jstar, k);
							*/				
							BasePair next = findNextCotranscriptionalBasePairMultipleHelix(canPair2, pairedProbs, diffs, Constants.IterationCutOff, sigma, jstar);
							if(next == null)
							{
								nextBasePair = getMaxBasePair(canPair2, pairedProbs, structure);
							}
							else
							{
								nextBasePair = new PPOutput(next.i, next.j, 0, new PointRes(0.0), new PointRes(0.0));								
							}
						}
						else
						{
							nextBasePair = getMaxBasePair(canPair2, pairedProbs, structure);
						}
						 
						
						
						
						int leftIdx = nextBasePair.getLeftIdx();
						int rightIdx = nextBasePair.getRightIdx();
						
						for(int i = 0 ; i < canPair2.length ; i++)
						{
							for(int j = i + 1 ; j < canPair2.length ; j++)
							{
								if(canPair2[i][j]  && structure[i] < 0 && structure[j] < 0)
								{
									if((leftIdx == -1 && rightIdx == -1) || pairedProbs[i][j].compareTo( pairedProbs[leftIdx][rightIdx]) > 0)
									{
										leftIdx = i;
										rightIdx = j;
									}
								}
							}
						}
						System.out.println("fast "+c+"\t"+leftIdx+"\t"+rightIdx);
						
						if ((leftIdx<0)||(rightIdx<0)) {
							ppProbs = new PPOutput(-1,-1,0,PointRes.ZERO,PointRes.ZERO);
						} else {
							boolean[][] incomp = new IncompatiblePairsFinder().find(canPair2, leftIdx, rightIdx);
							PointRes rprob = new IncompatiblePairsFinder().calculateComp(incomp,leftIdx,rightIdx,pairedProbs);
							PointRes[][] diffs2 = ppCalc.getDiffs(pairedProbs, unpairedProbs, canPair2);
							PointRes diff = diffs2[leftIdx][rightIdx];
						//	Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx,diffs,canPair2);
							Helix helix = new HelicesMaker().makeHelix(leftIdx,rightIdx, pairedProbs,diffs,canPair2);
							PPOutput ppProbs2 = new PPOutput(helix.getLeftIdx(), helix.getRightIdx(), helix.getHelixLength(), diff, rprob);
		
							//if(ppProbs2.getDiff().compareTo(new PointRes(Math.max(delta2,weight)))>0)
							
							// form the helix
							if(ppProbs2.getDiff().compareTo(new PointRes(Math.max(currentDelta*delta2,weight)))>0)
							{
								canPair = canPair2;
								ppProbs = ppProbs2;
								int [] structure2 = new int[structure.length];
								for(int i = 0 ; i < structure.length ; i++)
								{
									structure2[i] = structure[i]+1;
								}
								//System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(structure2));
								structure = new StructureUtils().makeNewStructure(structure, ppProbs2);
		
								continue;
							}
						}
						break;				
					}
				}
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			dumpStructure(structure);
		}


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		
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
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName()+metadata, structure);
	}
	
	
	public static BasePair findNextMaxPijBasePair (boolean [][] canPair, PointRes [][] pairedProbs, PointRes [][] diffs, PointRes deltaGreaterThan)
	{
		BasePair bp = null;
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = i + 1 ; j < canPair.length ; j++)
			{
				if(canPair[i][j] && diffs[i][j].compareTo(deltaGreaterThan) > 0)
				{
					if(bp == null || pairedProbs[i][j].compareTo(bp.probability) > 0)
					{
						bp = new BasePair(i, j, pairedProbs[i][j], diffs[i][j]);
					}
				}
			}
		}
		return bp;
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
	

	
	public static BasePair findNextCotranscriptionalBasePair(ArrayList<BasePair> possibleBasePairs, int jstar, double sigmaFactor)
	{
		PointRes sigma = new PointRes(sigmaFactor);
		
		BasePair minPair = null;
		PointRes min = null;
		int angelaOffset = 0;
		
		for(BasePair basePair : possibleBasePairs)
		{
			PointRes j =new PointRes(Math.max(basePair.i, basePair.j));
			//PointRes cotranscriptional = k.multiply(j);
			PointRes cotranscriptional = sigma.multiply(PointRes.max(PointRes.ZERO, j.subtract(new PointRes(jstar-angelaOffset))));

			if(basePair.probability.compareTo(PointRes.ZERO) == 0)
			{
				continue;
			}
			PointRes val = PointRes.ONE.divide(basePair.probability).add(cotranscriptional);
			//PointRes val = PointRes.ONE.divide(basePair.probability).add(k.multiply(cotranscriptional));
			if(min == null || val.compareTo(min) < 0)
			{
				minPair = basePair;
				min = val;
			}
		}
		
		return minPair;
	}
	
	public static BasePair findNextCotranscriptionalBasePairMultipleHelix(boolean [][] canPair, PointRes [][] pairedProbs, PointRes [][] diffs, PointRes deltaGreaterThan, double kfactor, int jstar)
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
}
