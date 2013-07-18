package uk.ac.ox.osscb.analysis;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;


public class CompareReversal {
	
	public static void main(String [] args) throws IOException{
		//open normal and reverse, through all in directory
		//reverse reversed set
		//count diff per each
		//use * on third line when printing out to represent different columns
		if(args.length == 2){
			// normal reversal
			String normalName = args[0];
			String reversalName = args[1]; 
			
			int diffCountTotal = 0; 
			
			File dataDir = new File("datasets/"); //deal with later
			String outputDirString = "output/";
			File outputNormalDir = new File(outputDirString + normalName);
			File outputReversalDir = new File(outputDirString + reversalName);
			String compareDirString = "Comparision/";
			File compareDir = new File(compareDirString + normalName + "_" + reversalName +"_comparision/");
			compareDir.mkdirs();
			
			BufferedWriter writerC = new BufferedWriter(new FileWriter(compareDir.getPath() + "/" +normalName + 
					"_" + reversalName + ".txt"));
			
			File[] outputReversalFiles = outputReversalDir.listFiles(new FilenameFilter() { 
   	         	public boolean accept(File dir, String filename)
   	         	{ return filename.endsWith(".dbn"); }
			} );
			
			for(File normalFile : outputNormalDir.listFiles()){
				
				if(normalFile.getName().endsWith(".dbn")){
					
					for (File revFile : outputReversalFiles){
						String basename = normalFile.getName().replace(".dat.fas.noevol.dbn", ""); 
						System.out.println(basename);
						if(revFile.getName().contains(basename+"_reversed")){
							System.out.println(normalFile.getName());
							String normStructure = findStructure(normalFile);
							String reversedStructure = reverse(findStructure(revFile));
							int count = 0; 
							if(normStructure.length() == reversedStructure.length()){
								for(int i = 0; i < normStructure.length(); i++){
									if(normStructure.charAt(i) != reversedStructure.charAt(i)){
										count++; 
										diffCountTotal++; 
									}
								}
							}
								
							writerC.write(basename + ":\r\n" + normStructure + "\r\n" + reversedStructure 
									+ "\r\ndiffCount = "+ count +"\r\n\r\n");
							break;
						}
					}
					
					
				}
			}
			writerC.write("Total difference Count = " + diffCountTotal + "\r\n");
	        writerC.close();
		}
	}
	
	public static String findStructure(File file) throws IOException{
		BufferedReader buffer = new BufferedReader(new FileReader(file));
		String textline = null;
		while ((textline = buffer.readLine()) != null) {
			if (!(textline.startsWith(">")) && (textline.contains("(") || textline.contains("."))) {
				buffer.close();
				return textline; 
			} 
	    }
		buffer.close();
	    return null; 
	}

	public static String reverse(String textline){
		StringBuffer reverseBuffer = new StringBuffer(textline);
    	reverseBuffer.reverse();
    	//"a" and "b" are intermediaries for reversing brackets
    	String tmp = reverseBuffer.toString().replace("(", "a").replace(")", "b");
    	return tmp.replace("a", ")").replace("b", "(");
	}
}
