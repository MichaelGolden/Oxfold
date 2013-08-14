package uk.ac.ox.osscb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputHelixInternalResult;
import uk.ac.ox.osscb.KineticFoldPppCalculatorBase.PPOutputInternalResult2;
import uk.ac.ox.osscb.analysis.IO;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.StructureData;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;
import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.GrammarParser;
import uk.ac.ox.osscb.grammar.RuleType;
import uk.ac.ox.osscb.inoutside.CoFoldInsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.CoFoldPosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.IOsideCalculator;
import uk.ac.ox.osscb.inoutside.InsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PPOutput;
import uk.ac.ox.osscb.inoutside.PPOutputHelix;
import uk.ac.ox.osscb.inoutside.ParallelInsideOutsideCalculator;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilities;
import uk.ac.ox.osscb.inoutside.PosteriorProbabilitiesCalculator;
import uk.ac.ox.osscb.inoutside.ValidatingIOSideCalculator;
import uk.ac.ox.osscb.parser.AlignmentParser;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;
import uk.ac.ox.osscb.phylo.PPfoldPhylogeneticCalculation;

public class CoFoldAnalogue2 {
	private final Logger log = LoggerFactory.getLogger(CoFoldAnalogue2.class);
	
	static int dummy = 0;
	
	public static boolean [] getGappyColumns(String [] align, double moreThanPercGaps)
	{
		double [] count = new double[align[0].length()];
		boolean [] gappy = new boolean[count.length];
		for(int i = 0 ; i < align.length ; i++)
		{
			for(int j = 0 ; j < align[0].length() ; j++)
			{
				if(align[i].charAt(j) == '-' || align[i].charAt(j) == 'N')
				{
					count[j]++;
				}
			}
		}
		
		for(int j = 0 ; j < count.length ; j++)
		{
			count[j] /= ((double)align.length);
			if(count[j] >= moreThanPercGaps)
			{
				gappy[j] = true;
			}
		}
		
		return gappy;
	}
	
	public static String [] deleteColumns(String [] align, boolean [] delete)
	{
		String [] ret = new String[align.length];
		
		double [] count = new double[delete.length];
		for(int i = 0 ; i < align.length ; i++)
		{
			ret[i] = "";
			for(int j = 0 ; j < count.length ; j++)
			{
				if(!delete[j])
				{
					ret[i] += align[i].charAt(j);
				}
			}
				
		}
		
		return ret;
	}
	
	public static boolean [][] deleteColumns(boolean [][] matrix, boolean [] delete)
	{
		ArrayList<boolean[]> list = new ArrayList<boolean[]>();
		
		for(int i = 0 ; i < matrix.length ; i++)
		{
			if(!delete[i])
			{
				list.add(matrix[i]);
			}
		}
		
		int cols = list.size();
		boolean [][] ret = new boolean[cols][cols];
		for(int i = 0 ; i < list.size() ; i++)
		{
			boolean [] oldrow = list.get(i); 
			boolean [] newrow = new boolean[cols];
			int k = 0;
			for(int j = 0 ; j < oldrow.length ; j++)
			{
				if(!delete[j])
				{
					newrow[k] = oldrow[j];
					k++;
				}
			}			
			ret[i] = newrow;
		}
		
		return ret;			
	}
	
	public static void printMatrix(boolean [][] matrix)
	{
		for(int i = 0 ; i < matrix.length ; i++)
		{
			for(int j = 0 ; j < matrix[0].length ; j++)
			{
				System.out.print(matrix[i][j] ? 1 : 0);
			}

			System.out.println();
		}
	}
	
	public static int [] reinsertDeleted2(int [] structure, boolean [] deleted)
	{
		String dbs = RNAFoldingTools.getDotBracketStringFromPairedSites(structure);
		String dbs2="";
		int [] ret = new int[deleted.length];
		int a = 0;
		for(int i = 0 ; i < ret.length ; i++)
		{	
			if(!deleted[i])
			{
				dbs2 += dbs.charAt(a);
				a++;
			}
			else
			{
				dbs2 += ".";
			}
		}
		
		return RNAFoldingTools.getPairedSitesFromDotBracketString(dbs2);
	}
	
	public static int [] reinsertDeleted(int [] structure, boolean [] deleted, int unpairedBase)
	{
		int a = 0;
		int b = 0;
		int [] s = new int[deleted.length];
		for (int i = 0; i < deleted.length; i++) {
			s[i] = unpairedBase;
			if (!deleted[i]) {
				if (structure[a] != unpairedBase) {
					s[i] = structure[a] + b;
				}
				a++;
			} else {
				b++;
			}
		}
		return s;
	}
	
	
	
	/*
	public static int [] reinsertDeleted(int [] structure, boolean [] deleted, int unpairedbase)
	{
		int [] ret = new int[deleted.length];
		int a = 0;
		for(int i = 0 ; i < ret.length ; i++)
		{
			ret[i] = unpairedbase;
			if(!deleted[i])
			{
				ret[i] = structure[a];
				a++;
			}
		}
		return ret;
	}
	*/
	
	public static boolean [][] reinsertDeleted(boolean [][] matrix, boolean [] deleted)
	{
		boolean [][] ret = new boolean[deleted.length][deleted.length];

		int a = 0;
		for(int i = 0 ; i < deleted.length ; i++)
		{
			int b = 0;
			for(int j = 0 ; j < deleted.length ; j++)
			{						
				if(!deleted[j])
				{
					if(!deleted[i])
					{
						ret[i][j] = matrix[a][b];
					}
					b++;
				}
			}
			
			if(!deleted[i])
			{
				a++;
			}
		}
		
		return ret;
	}
	
	public static PointRes [][] reinsertDeleted(PointRes [][] matrix, boolean [] deleted, PointRes emptyVal)
	{
		PointRes [][] ret = new PointRes[deleted.length][deleted.length];

		int a = 0;
		for(int i = 0 ; i < deleted.length ; i++)
		{
			int b = 0;
			for(int j = 0 ; j < deleted.length ; j++)
			{		
				ret[i][j] = emptyVal;
				if(!deleted[j])
				{
					if(!deleted[i])
					{
						ret[i][j] = matrix[a][b];
					}
					b++;
				}
			}
			
			if(!deleted[i])
			{
				a++;
			}
		}
		
		return ret;
	}
	
	public static void main(String [] args)
	{
		boolean [][] matrix = {{true, true, true,false},{true, false,false,true},{true,true,false, false},{true,false,true,true}};
		printMatrix(matrix);
		boolean [] delete = {false,true,false,false};
		boolean [][] ret = deleteColumns(matrix,delete);		
		System.out.println();
		printMatrix(ret);
		boolean [][] ret2 = reinsertDeleted(ret, delete);
		printMatrix(ret2);
		
		
		int [] structure = {1,2,3};
		int [] structure2 = reinsertDeleted(structure, delete, -1);
		for(int i = 0 ; i < structure2.length ; i++)
		{
			System.out.println(structure2[i]);
		}
	}
	
	
	public static NucleotideProbsPrecise cotranscriptionalWeight(NucleotideProbsPrecise nucleotideProbs, double alpha, int numSeq)
	{
		PointRes [][] probs = nucleotideProbs.getMtx();
		PointRes [] unpairing = nucleotideProbs.getVector();
		PointRes [][] sum = new PointRes[probs.length][probs.length];
		
		int maxDistance = 100000000;
		
		for(int i = 0 ; i < unpairing.length ; i++)
		{
			probs[i][i] = unpairing[i];
		}
		
		for(int i = 0 ; i < sum.length ; i++)
		{			
			for(int j = 0 ; j < sum.length ; j++)
			{			
				if(probs[i][j] != null)
				{
					probs[j][i] = probs[i][j];
				}
			}
		}
		
		PointRes minProbability = PointRes.valueOf(0.0);
		for(int i = 0 ; i < sum.length ; i++)
		{			
			sum [i][0] = PointRes.ZERO;
			for(int j = 0 ; j <= i ; j++)
			{
				sum[i][j] = PointRes.ZERO;
			}
			for(int j = 1  ; j < sum.length ; j++)
			//for(int j = i+1 ; j < sum.length ; j++)
			{				
				if(i != j -1)
				{					
					//if((i >= j || Math.abs(i-j) < maxDistance) && (all || probs[i][j].compareTo(PointRes.valueOf(0.1)) > 0))
					if(probs[i][j-1].compareTo(minProbability) > 0 && Math.abs(i-j) < maxDistance)
					{
						sum[i][j] = sum[i][j-1].add(probs[i][j-1].doubleValue());
					//	sum[i][j] = sum[i][j-1].add(Math.pow(probs[i][j-1].doubleValue(), 1/((double)numSeq)));
					}
					else
					{
						sum[i][j] = sum[i][j-1];
					}
				}
				else
				{
					sum[i][j] = sum[i][j-1];
				}
			}
		}
		
		
		
		double max = Math.max(0, 1 - alpha);
		double beta = 50;
		NucleotideProbsPrecise nucProbs = new NucleotideProbsPrecise(probs.length);
		PointRes [][] ret = new PointRes[probs.length][probs.length];		
		for(int i = 0 ; i < ret.length ; i++)
		{
			PointRes sumA = PointRes.ZERO;
			PointRes sumB = PointRes.ZERO;
			ret[i][i] = probs[i][i];
			for(int j = i+1 ; j < ret.length ; j++)
			{
				//if(i != j && (i >= j || Math.abs(i-j) < maxDistance) && (all || probs[i][j].compareTo(PointRes.valueOf(0.1)) > 0))
				if(probs[i][j].compareTo(minProbability) > 0  && Math.abs(i-j) < maxDistance)
				{
					double normalizedSum = (sum[i][j].doubleValue() / (double)(sum.length))*beta;
					ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(normalizedSum)).multiply(probs[i][j]);
				//	ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(sum[i][j].doubleValue())).multiply(probs[i][j]);		
					//ret[i][j] = PointRes.valueOf(1 + max*Math.tanh(beta*((sum[i][ret.length-1].doubleValue() - sum[i][j].doubleValue()))/ (double)(sum.length))).multiply(ret[i][j]);
					ret[j][i] = ret[i][j];		
					
					//ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(sum[i][j].doubleValue())).multiply(probs[i][j]);					
					//ret[i][j] = PointRes.valueOf(1 + max*Math.tanh(sum[i][ret.length-1].doubleValue() - sum[i][j].doubleValue())).multiply(ret[i][j] );
		
					
					sumA = sumA.add(ret[i][j]);
					sumB = sumB.add(probs[i][j]);
					
					//ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(sum[i][j].doubleValue())).multiply(PointRes.ONE);
					//System.out.println("XX"+Math.tanh(sum[i][j].doubleValue()));
					//System.out.println("XY"+Math.tanh(sum[i][ret.length-1].doubleValue() - sum[i][j].doubleValue()));
					
					//System.out.println(i+"\t"+j+"\t"+ret[i][j].doubleValue());
					//double normalizedSum = (sum[i][j].doubleValue() / (double)(sum.length))*beta;
					//ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(normalizedSum)).multiply(probs[i][j]);
					//ret[i][j] = PointRes.valueOf(1 + max*Math.tanh(beta*((sum[i][ret.length-1].doubleValue() - sum[i][j].doubleValue()))/ (double)(sum.length))).multiply(ret[i][j]);
//					ret[i][j] = PointRes.valueOf(1 - max*Math.tanh(sum[i][j].doubleValue())).multiply(probs[i][j]);
		
					
				}
				else
				{	
					ret[i][j] = probs[i][j];
					ret[j][i] = ret[i][j];
				}
			}


			if(sumA.compareTo(PointRes.ZERO) != 0)
			{
				PointRes normalisation = sumB.divide(sumA);
				for(int j = i+1 ; j < ret.length ; j++)
				{
					//if(i != j  && (i >= j || Math.abs(i-j) < maxDistance) && (all || probs[i][j].compareTo(PointRes.valueOf(0.1)) > 0))
					if(probs[i][j].compareTo(minProbability) > 0  && Math.abs(i-j) < maxDistance)
					{
						ret[i][j] = ret[i][j].multiply(normalisation);
						ret[j][i] = ret[i][j];		
					}
				}
			}
		}
		
		for(int i = 0 ; i < ret.length ; i++)
		{
			nucProbs.setUnpairingProbability(i, nucleotideProbs.getUnpairingProbability(i));
			for(int j = 0  ; j < ret.length ; j++)
			{
				nucProbs.setPairingProbability(i, j, ret[i][j]);
			}
		}
		
		return nucProbs;
	}
	
	public void foldEvolutionary(String alignmentFile, String grammarFile, String paramsFile, String treeFile, double alpha, double tau){
		Util.assertCanReadFile(alignmentFile);
		Util.assertCanReadFile(grammarFile);
		Util.assertCanReadFile(paramsFile);
		Util.assertCanReadFile(treeFile);
		/*if(weight < 0){
			throw new IllegalArgumentException(String.format("Weight must be non-negative. Input: %f", weight));
		}*/

		Grammar grammar = GrammarParser.parse(grammarFile);
		
		EvolutionaryParameters parameters = ParameterParserEvolutionary.parse(paramsFile);
		
		AlignmentParser alignParse = new DefaultAlignmentParser();
		
		String[] align = alignParse.parseEvolutionary(alignmentFile,parameters.getSAlphabet());
		
		
		/*
		boolean [] delete = CoFoldAnalogue.getGappyColumns(align, 0.25);
		for(int j = 0 ; j < align.length; j++)
		{
			for(int i = 0 ; i < align[j].length() ; i++)
			{
				System.out.print(align[j].charAt(i));
			}
			System.out.println();
			for(int i = 0 ; i < delete.length ; i++)
			{
				System.out.print(delete[i] ? 1 : 0);
			}
			System.out.println();
		}
		System.out.println();
		
		//boolean [] delete = new boolean[align[0].length()];
		//boolean [] delete2 = CoFoldAnalogue.getGappyColumns(align, 0.75);
		*/
		//boolean [] delete = new boolean[align[0].length()];
		boolean [] delete = CoFoldAnalogue2.getGappyColumns(align, 0.25);
		align = CoFoldAnalogue2.deleteColumns(align, delete);
		
		/*
		for(int j = 0 ; j < align.length; j++)
		{
			for(int i = 0 ; i < align[j].length() ; i++)
			{
				System.out.print(align[j].charAt(i));
			}
			System.out.println();
			for(int i = 0 ; i < delete.length ; i++)
			{
				System.out.print(delete[i] ? 1 : 0);
			}
			System.out.println();
		}
		System.out.println();*/
		
				
		EvolutionaryTree tree = new EvolutionaryTreeParser().parse(treeFile);
		
		NucleotideProbsPrecise alignmentProbs = new NucleotideBasePairingProbsCalculator().calculate(align, parameters);
		NucleotideProbsPrecise alignmentProbsEvol = new EvolProbabilitiesCalculator().getEvolutionaryProbs(tree, parameters, align);

		
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();
		IO.loadFastaSequences(new File(alignmentFile), sequences, sequenceNames);
		alignmentProbs.writeEvolutionaryProbs(new File(alignmentFile+".nonevolutionary"));

		//alignmentProbsEvol = alignmentProbs;
	
		//NucleotideProbsPrecise cotranscriptional = cotranscriptionalWeight(alignmentProbsEvol);


		alignmentProbsEvol = EvolProbabilitiesCalculator.calculateEvolutionaryProbsPPfold(align, sequenceNames, new File(treeFile));
		alignmentProbsEvol.writeEvolutionaryProbs(new File(alignmentFile+".evolutionary"));
		
		NucleotideProbsPrecise cotranscriptional = cotranscriptionalWeight(alignmentProbsEvol, 0.1, 5);
		alignmentProbsEvol = cotranscriptional;
		//NucleotideProbsPrecise cotranscriptional = alignmentProbsEvol;
		cotranscriptional.writeEvolutionaryProbs(new File(alignmentFile+".cotranscriptional"));
		System.out.println("alpha="+alpha);
		
		if(log.isDebugEnabled()){
			log.debug(String.format("Pairing probs (NON-evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (NON-evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
			
			log.debug(String.format("Pairing probs (Evol):\r\n%s", Util.print2DArray(alignmentProbs.getMtx())));
			log.debug(String.format("Unpairing probs (Evol):\r\n%s", Util.dump1DArray(alignmentProbs.getVector(), 2, 4)));
		}


		//CoFoldInsideOutsideCalculator ioCalc = new CoFoldInsideOutsideCalculator(grammar);
		ParallelInsideOutsideCalculator parallelCalc = new ParallelInsideOutsideCalculator(grammar);
		
		
		//IOsideCalculator ioCalc = new ValidatingIOSideCalculator(new CoFoldInsideOutsideCalculator(grammar), grammar);
		//alpha = 0.5;
		//tau = 640;
		//CoFoldPppCalculator cofoldCalc = new CoFoldPppCalculator(alpha, tau,grammar,ioCalc);
		/*KineticFoldPppCalculator kFPppCalc = weight > 0 ? new KineticFoldPppCalculatorWithWeight(weight, grammar, ioCalc)
			// 0 == weight, negative was rejected at the very beginning
			:new KineticFoldPppCalculatorWeightLess(grammar, ioCalc);*/
		
		//ioCalc.outside(insideProbs, pairingProbs, distances, alpha, tau, structure, canPair)
		
		boolean evol = true;
		boolean exitBecauseOfDiff = false;
	
		boolean [][] canPair = new boolean [delete.length][delete.length];
		for(int i = 0 ; i < canPair.length ; i++)
		{
			for(int j = 0 ; j < canPair.length ; j++)
			{
				if(Math.abs(i-j) > 3)
				{
					canPair[i][j] = true;
				}
			}	
		}
		canPair = CoFoldAnalogue2.deleteColumns(canPair, delete);
				
		// by default is initialised with zeros automatically
		int[] structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		//InsideOutsideProbabilities insideProbs = ioCalc.inside(evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		//InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		long startIOSerial = System.nanoTime();
		//InsideOutsideProbabilities insideProbs = ioCalc.inside(evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		//InsideOutsideProbabilities outsideProbs = ioCalc.outside(insideProbs, evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		long endIOSerial = System.nanoTime();
		System.out.println("IO serial: "+((endIOSerial - startIOSerial) / 1e9));
		
		structure = new int[align[0].length()];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
		long startIOParallel = System.nanoTime();
		InsideOutsideProbabilities insideProbs = parallelCalc.insideParallelCoFold(evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		InsideOutsideProbabilities outsideProbs = parallelCalc.outsideParallelCoFold(insideProbs, evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		//insideProbs = parallelCalc.insideParallelCoFold(evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		//outsideProbs = parallelCalc.outsideParallelCoFold(insideProbs, evol ? alignmentProbsEvol : alignmentProbs, alpha, tau, structure, canPair);
		
		long endIOParellel = System.nanoTime();
		//System.out.println("IO parallel: "+((endIOParellel - startIOParallel) / 1e9));
		double speedup =((endIOSerial - startIOSerial) / 1e9)/((endIOParellel - startIOParallel) / 1e9);
		//System.out.println("Speedup "+speedup);
				
		CoFoldPosteriorProbabilitiesCalculator ppCalc = new CoFoldPosteriorProbabilitiesCalculator(grammar);
		int [][] distances = null;
		PosteriorProbabilities currentPostProbs = ppCalc.calculate(insideProbs, outsideProbs, evol ? alignmentProbsEvol : alignmentProbs, distances, alpha, tau, structure, canPair);
		//MEACalculator meaCalculator = new MEACalculator();
		int [] decoded_structure = RNAFoldingTools.getPosteriorDecodingConsensusStructure(currentPostProbs.getBasePairProbs());
		currentPostProbs.pairedProbs = CoFoldAnalogue2.reinsertDeleted(currentPostProbs.pairedProbs, delete, PointRes.ZERO);
		
		
		//structure = structure2;
		decoded_structure = CoFoldAnalogue2.reinsertDeleted(decoded_structure, delete, 0);
		String dbs = RNAFoldingTools.getDotBracketStringFromPairedSites(decoded_structure);
		structure = new int[decoded_structure.length];
		for(int i = 0 ; i < structure.length ; i++)
		{
			structure[i] = decoded_structure[i]-1;
		}
		//structure = meaCalculator.calculate(currentPostProbs);
		
		/*
		for (ProductionRule rule1: grammar.getRules(RuleType.RULE1)) {
			try {
				outsideProbs.writeTable(new File("oprobs_"+rule1.getLeft()+"_"+dummy+".txt"), rule1.getLeft());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		
		for (ProductionRule rule2: grammar.getRules(RuleType.RULE2)) {
			try {
				outsideProbs.writeTable(new File("oprobs_"+rule2.getLeft()+"_"+dummy+".txt"), rule2.getLeft());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		for (ProductionRule rule3: grammar.getRules(RuleType.RULE3)) {
			try {
				outsideProbs.writeTable(new File("oprobs_"+rule3.getLeft()+"_"+dummy+".txt"), rule3.getLeft());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		dummy++;*/
		
		//
		//System.out.println(RNAFoldingTools.);
		
		//boolean[][] canPair = new PossiblePairFinder().canPair(structure);
		//InsideOutsideProbabilities insideProbs = ioCalc.insideE(alignmentProbs, structure, canPair);
		//InsideOutsideProbabilities outsideProbs = ioCalc.outsideE(insideProbs, alignmentProbs, structure, canPair);
		//PosteriorProbabilitiesCalculator ppCalc = new PosteriorProbabilitiesCalculator(grammar);
		//PosteriorProbabilities currentPostProbs = ppCalc.calculateE(insideProbs, outsideProbs, alignmentProbs, structure, canPair);
		/*
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			canPair = new PossiblePairFinder().canPair(structure);			
					
			PPOutputInternalResult2 postProbs = evol ? 
					cofoldCalc.calculatePpOutputInternalResult2(alignmentProbsEvol, structure, currentPostProbs) 
						: cofoldCalc.calculatePpOutputInternalResult2(alignmentProbs, structure, currentPostProbs);
			
			PPOutput ppProbs = postProbs.getPpProbs();
			currentPostProbs = postProbs.getPosteriorProbs();

			//if (ppProbs.getDiff().signum()>0) {
			if ((!evol)&&(ppProbs.getDiff().compareTo(Constants.EndNonEvolutionaryFold)<0)) {
				evol = true;
				continue;
			}
			if (ppProbs.getDiff().compareTo(Constants.IterationCutOff)>0) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs);
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			// iterationsGenerator.generate(structure);
			// OutputGenerator outputGenerator = new LoggingOutputGenerator();
			// outputGenerator.generate(structure);
			dumpStructure(structure);
		}*/


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);
		outputGenerator.generateFinal(structure);
		
		try {
			currentPostProbs.savePosteriorProbabilities(new File(alignmentFile+".evol.bp"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		writeDotBracketFile(new File(alignmentFile+".evol.dbn"),new File(alignmentFile).getName(), structure);
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
				, ppProbs.getComp()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}
}
