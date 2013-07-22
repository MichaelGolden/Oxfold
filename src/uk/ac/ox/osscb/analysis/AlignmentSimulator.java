package uk.ac.ox.osscb.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class AlignmentSimulator {
	
	public static String alphabet = "ACGT";
	public static double [] single = {0.364097, 0.273013, 0.211881, 0.151009};	


	public static double [][] paired = {{0.001167, 0.001806, 0.001058, 0.177977}
						,{0.001806 , 0.000391, 0.266974, 0.000763}
						,{0.001058 , 0.266974, 0.000406, 0.049043}
						,{0.177977 , 0.000763, 0.049043, 0.002793}};
	
	public static char getUnpaired(double beta)
	{		
		double sum = 0;
		for(int i = 0 ; i < single.length ; i++)
		{
			if(beta >= sum && beta < sum + single[i])
			{
				return alphabet.charAt(i);
			}
			
			sum += single[i];
		}
		
		return 'N';
	}
	
	public static String getPaired(double beta)
	{		
		double sum = 0;
		for(int i = 0 ; i < paired.length ; i++)
		{	
			for(int j = 0 ; j < paired.length ; j++)
			{
				if(beta >= sum && beta < sum + paired[i][j])
				{
					return ""+alphabet.charAt(i)+alphabet.charAt(j);
				}
				
				sum += paired[i][j];
			}
		}
		
		return "NN";
	}
	
	public static String simulate(Random random, String sequence, int [] pairedSites, double expectedMutations)
	{
		
		StringBuffer sb = new StringBuffer(sequence);
		for(int i = 0 ; i < sequence.length() ; i++)
		{
			double alpha = pairedSites[i] == 0 ? expectedMutations :  expectedMutations/2;
			if(random.nextDouble() < alpha)
			{
				double beta = random.nextDouble();
				if(pairedSites[i] == 0)
				{
					sb.setCharAt(i, getUnpaired(beta));
				}
				else
				{
					String pair = getPaired(beta);
					sb.setCharAt(i, pair.charAt(0));
					sb.setCharAt(pairedSites[i]-1, pair.charAt(1));
				}
			}
		}	
		
		return sb.toString().replace("T", "U");		
	}
	
	public static void simulate(ArrayList<String> sequences,  ArrayList<String> sequenceNames, Random random, String sequence, int [] pairedSites, double expectedMutations, int n, boolean includeFirst)
	{
		ArrayList<String> retSequences = simulate(random, sequence, pairedSites, expectedMutations, n, includeFirst);

		for(int i = 0 ; i <retSequences.size() ; i++)
		{
			sequences.add(retSequences.get(i));
			sequenceNames.add("seq"+i);
		}
	
	}
	
	public static ArrayList<String> simulate(Random random, String sequence, int [] pairedSites, double expectedMutations, int n, boolean includeFirst)
	{
		ArrayList<String> sequences = new ArrayList<String>();
		if(includeFirst)
		{
			sequences.add(sequence);
		}
		while(sequences.size() < n)
		{
			sequences.add(simulate(random, sequence, pairedSites, expectedMutations));
		}
		return sequences;
	}
	
	public static void main(String [] args)
	{
		
		try {
			StructureData data = StructureData.readExperimentalStructureData(new File("datasets/PDB_00488_23s_rRNA_HALOARCULA_MARISMORTUI.txt_0_527.dat"));

			Random random = new Random(3494163181355294561L);
			String sequence = data.sequences.get(0);
			int [] pairedSites = data.pairedSites;
			ArrayList<String> sequences = simulate(random, sequence, pairedSites, 0.2, 5, true);

			System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites));
			for(String s : sequences)
			{
				System.out.println(s);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
