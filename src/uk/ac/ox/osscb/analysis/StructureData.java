package uk.ac.ox.osscb.analysis;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class StructureData implements Comparable<StructureData> {
	
	public File file;
	public String title="";
	
	public int [] pairedSites; // structure
	
	public File basePairProbFile;
	public double [][] basePairProb;
	
	public ArrayList<String> sequences = new ArrayList<String>();
	public ArrayList<String> sequenceNames = new ArrayList<String>();
	
	double time;
	
	public StructureData()
	{
		
	}
	
	public StructureData(int [] pairedSites)
	{
		this.pairedSites = pairedSites;
	}
	
	/*
	 * Read a file containing an experimentally determined structure and an alignment.
	 * 1st line contains structure in dot bracket format (using angle brackets).
	 * Remaining lines contain alignment in FASTA format.
	 */
	public static StructureData readExperimentalStructureData(File experimentalFile) throws IOException
	{
		StructureData data = new StructureData();
		data.file = experimentalFile;
		BufferedReader buffer = new BufferedReader(new FileReader(experimentalFile));
		data.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(buffer.readLine().replace('<', '(').replace('>', ')'));
		String textline = null;
		String sequence = "";
        
        while ((textline = buffer.readLine()) != null) {
            if (textline.startsWith(">")) {
                data.sequenceNames.add(textline.substring(1));
                if (!sequence.equals("")) {
                	data.sequences.add(sequence.replace(".", "-").toUpperCase());
                    sequence = "";
                }
            } else {
                sequence += textline.trim();
            }

        }
        buffer.close();
        if (!sequence.equals("")) {
        	data.sequences.add(sequence.replace(".", "-").toUpperCase());
        }
        /*
        ArrayList<String> sequences = new ArrayList<String>();
        ArrayList<String> sequenceNames = new ArrayList<String>();
        String seq = data.sequences.get(0);
        sequences.add(data.sequences.get(0).replaceAll("-", ""));
        sequenceNames.add(data.sequenceNames.get(0));
        data.sequences = sequences;
        data.sequenceNames = sequenceNames;
        

        int [] s = new int[sequences.get(0).length()];

        
    	int a = 0;
    	int b = 0;
        for(int i = 0 ; i < data.pairedSites.length ; i++)
        {        	
        	if(seq.charAt(i) != '-')
        	{
        		if(data.pairedSites[i] != 0)
        		{
        			s[a] = data.pairedSites[i]-b;
        		}
        		//System.out.println(a+"\t"+s[a]);
        		a++;
        	}
        	else
        	{
        		b++;
        	}
        }
        
        System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(data.pairedSites));
        System.out.println(RNAFoldingTools.getDotBracketStringFromPairedSites(s));
        
        data.pairedSites = s;
        */
        
        return data;
	}
	
	public String toString()
	{
		String ret = "";
		ret += RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites)+"\n";
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			ret += ">"+sequenceNames.get(i)+"\n";
			ret += sequences.get(i)+"\n";
		}
		return ret;		
	}

	@Override
	public int compareTo(StructureData o) { // sort by paired sites lengths
		return this.pairedSites.length - o.pairedSites.length;
	}
	
	public static double [] getNucleotidePairedProbs(double [][] basePairProb)
	{
		double [] probs = new double[basePairProb.length];
		for(int i = 0 ; i < basePairProb.length ; i++)
		{
			for(int j = 0 ; j < basePairProb.length ; j++)
			{
				probs[i] += basePairProb[i][j];
			}
		}
		return probs;
	}
	
	public static double [][] getBasePairProb(File bpFile) throws IOException
	{
		double [][] basePairProb = null;
		
		BufferedReader buffer = new BufferedReader(new FileReader(bpFile));
		String textline = null;
		int i = 0;
		while((textline = buffer.readLine()) != null)
		{
			String [] split = textline.split("[,(\\s)](\\s)*");
			if(i == 0)
			{
				basePairProb = new double[split.length][split.length];
			}
			
			for(int j = 0 ; j < basePairProb.length ; j++)
			{
				basePairProb[i][j] = Double.parseDouble(split[j]);
			}
			
			i++;
		}
		
		for(int x = 0 ; x < basePairProb.length ; x++)
		{
			for(int y = 0 ; y < basePairProb.length ; y++)
			{
				if(basePairProb[x][y] != 0)
				{
					basePairProb[y][x] = basePairProb[x][y];
				}
			}
		}
		
		buffer.close();
		
		return basePairProb;
	}
	
}
