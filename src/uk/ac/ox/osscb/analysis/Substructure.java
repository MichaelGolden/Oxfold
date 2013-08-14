package uk.ac.ox.osscb.analysis;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Class that represents a substructure
 *
 * @author Michael Golden
 */
public class Substructure implements Serializable {

    public int[] pairedSites;
    public String name = "";
    /**
     * String representing the nucleotide sequence of this structure.
     */
    public String sequence = "";
    public int startPosition = 0;
    public int length;
    public int index;

    public Substructure(int length) {
        this.length = length;
        this.pairedSites = new int[length];
    }
    
    public Substructure(int startPosition, int [] pairedSites)
    {
        this.startPosition = startPosition;
        this.length = pairedSites.length;
        this.pairedSites = pairedSites;
    }

    public String getDotBracketString() {
        return RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites);
    }

    /**
     *
     * @return the start position of the structure in the parent genome (one-offset)
     */
    public int getStartPosition ()
    {
        return startPosition;
    }

    /**
     * @return the end position (inclusive) of the structure in the parent genome.
     */
    public int getEndPosition ()
    {
        return startPosition+length;
    }

    public int getLength ()
    {
        return length;
    }

    @Override
    public String toString ()
    {
        return name + " ["+getStartPosition()+"-"+getEndPosition()+" ("+length+")]\n"+(sequence == null ? "" : sequence+"\n")+getDotBracketString();
    }
    
    

    /*public static void main(String [] args)
    {
        try {
            Structure s = StructureParser.parseNaspCtFile(new File("D:/NASP/BFDV/BFDV_10Seq.ct"));
            for (int i = 0; i < s.length; i++)
            {
                System.out.println(s.pairedSites[0][i]+"\t"+s.pairedSites[1][i]);
            }
            System.out.println(s.getDotBracketString());
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }*/

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Substructure other = (Substructure) obj;
        if (this.startPosition != other.startPosition) {
            return false;
        }
        if (this.length != other.length) {
            return false;
        }
        return true;
    }
    
    public static void main(String [] args)
    {
    	//File rnastrandfile=new File("C:/Users/Michael/Dropbox/Oxfold II/datasets/PDB_00488_23s_rRNA_HALOARCULA MARISMORTUI.txt");
    	
    	File [] files = new File("C:/Users/Michael/Dropbox/Oxfold II/long_datasets/").listFiles();
    	for(File rnastrandfile : files)
    	{
    		if(!rnastrandfile.getName().endsWith(".dat"))
    		{
		    	Substructure structure = getRNAstrandStructure(rnastrandfile);
		    	//System.out.println(s);
		    	int [] pairedSites = structure.pairedSites;
		    	
		    	ArrayList<Substructure> substructures = enumerateAdjacentSubstructures(pairedSites,structure.sequence,0,1800,false);
		    	System.out.println("length="+pairedSites.length);
		    	for(int i = 0 ; i < substructures.size() ; i++)
		    	{
		    		System.out.println(substructures.get(i)+"\n");
		    		Substructure s = substructures.get(i);
		    		/*String name = rnastrandfile.getName()+"_"+s.startPosition+"_"+(s.startPosition+s.length)+".dat";
		    		saveDatFile(substructures.get(i),s.startPosition+"_"+(s.startPosition+s.length),new File(rnastrandfile.getParentFile().getAbsolutePath()+File.separatorChar+name));
		    	}*/}
		    	String name = rnastrandfile.getName()+"_full.dat";
		    	//Substructurenew Substructure(0, pairedSites)
		    	saveDatFile(structure,name,new File(rnastrandfile.getParentFile().getAbsolutePath()+File.separatorChar+name));
    		}
    	}
    	
    	

    }
    
    public static void saveDatFile(Substructure s, String name, File outFile)
    {
    	try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
			writer.write(RNAFoldingTools.getDotBracketStringFromPairedSites(s.pairedSites).replace('(', '<').replace(')', '>').substring(0, s.sequence.length()));
			writer.newLine();
			writer.write(">"+name);
			writer.newLine();
			writer.write(s.sequence);
			writer.newLine();
			writer.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
    
    public static Substructure getRNAstrandStructure(File file)
    {
    	try {
			BufferedReader buffer = new BufferedReader(new FileReader(file));
			String textline = null;

			String header="";
			String sequence="";
			String structure="";
			int phase = 0;
			while((textline = buffer.readLine()) != null)
			{
				if(phase < 3 && textline.startsWith("#"))
				{
					header += textline+"\n";
					phase = 1;
				}
				else
			    if(phase < 3 && textline.matches("^[ACGTU]+$"))
				{
			    	phase = 2;
			    	sequence += textline;
				}
			    else
			    if(phase >= 2)
			    {
			    	phase = 3;
			    	structure += textline;
			    }
			}
			buffer.close();
			/*
			System.out.println(header.trim());
			System.out.println("----");
			System.out.println(sequence.trim());
			System.out.println("----");
			System.out.println(structure.trim());
			System.out.println("----");*/
			
			Substructure s = new Substructure(sequence.length());
			s.sequence = sequence;
			s.pairedSites = RNAFoldingTools.getPairedSitesFromDotBracketString(structure);
			return s;
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	return null;
    }
    
    public static ArrayList<Substructure> enumerateAdjacentSubstructures(int[] pairedSitesIn, String sequence, int minLength, int maxLength, boolean circularize) {
        ArrayList<Substructure> structures = new ArrayList<Substructure>();


        int[] pairedSites = Arrays.copyOf(pairedSitesIn, pairedSitesIn.length);
        int genomeLength = pairedSites.length;
        for (int i = 0; i < pairedSites.length; i++) {
            int x = i;
            int y = pairedSites[i];

            if (y > 0 && y - x + 1 > 0) {
                int[] pairedSitesSub = new int[y - x];
                for (int j = 0; j < pairedSitesSub.length; j++) {
                    if (pairedSites[i + j] != 0) {
                        pairedSitesSub[j] = pairedSites[i + j] - i;
                    } else {
                        pairedSitesSub[j] = 0;
                    }
                }

                Substructure s = new Substructure(y - x);
                s.pairedSites = pairedSitesSub;
                s.startPosition = x;
                s.name = structures.size() + "";

                if (maxLength == 0 || s.length <= maxLength) {
                    i += s.length;
                    if (s.length >= minLength && x + s.length < genomeLength) {
                        System.out.println(structures.size() + "\t" + s.getDotBracketString());
                        structures.add(s);
                    }
                }
            }
        }
        
        for(Substructure s : structures)
        {
        	
        	while(s.startPosition > 0 && pairedSites[s.startPosition-1] == 0)
        	{
        		s.startPosition--;
            	s.length++;
        	}
        	
        	while(s.startPosition+s.length < pairedSites.length && pairedSites[s.startPosition+s.length] == 0)
        	{
        		s.length++;
        	}
        	
        	s.pairedSites = new int[s.length];
        	for(int i = 0 ; i < s.pairedSites.length ; i++)
        	{
        		if(pairedSites[s.startPosition+i] != 0)
        		{
        			s.pairedSites[i] = pairedSites[s.startPosition+i] - s.startPosition;        		
        		}
        	}
        	
        	if(sequence != null)
			{
        		//System.out.println(sequence.length()+"\t"+s.startPosition+"\t"+Math.min(s.startPosition+s.length, sequence.length()));
				s.sequence = sequence.substring(s.startPosition,Math.min(s.startPosition+s.length, sequence.length()));
			}
        }

        return structures;
    }
}
