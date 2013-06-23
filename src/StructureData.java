import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class StructureData {
	
	File file;
	
	int [] pairedSites; // structure
	
	ArrayList<String> sequences = new ArrayList<>();
	ArrayList<String> sequenceNames = new ArrayList<>();
	
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
                	data.sequences.add(sequence.toUpperCase());
                    sequence = "";
                }
            } else {
                sequence += textline.trim();
            }

        }
        buffer.close();
        if (!sequence.equals("")) {
        	data.sequences.add(sequence.toUpperCase());
        }
        
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
}
