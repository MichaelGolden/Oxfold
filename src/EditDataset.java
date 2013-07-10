import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class EditDataset {
	private static File dataDir = new File("datasets/");
	
	public static void main(String [] args) throws IOException{
		
		for(File experimentalFile : dataDir.listFiles()){
			reverseDataset(experimentalFile.getName());
		}
		
	}
	
	public static void reverseDataset(String filename) throws IOException{
		BufferedReader buffer = new BufferedReader(new FileReader(dataDir + "/" + filename));
		String textline = null;
		BufferedWriter writerRD = new BufferedWriter(new FileWriter(dataDir + "/"  + filename.replace(".dat", "") + "_reversed.dat"));
		
        while ((textline = buffer.readLine()) != null) {
            if (!(textline.startsWith(">"))) {
            	StringBuffer reverseBuffer = new StringBuffer(textline);
            	reverseBuffer.reverse();
            	//"a" and "b" are intermediaries for reversing brackets
            	String tmp = reverseBuffer.toString().replace("<", "a").replace(">", "b");
            	tmp = tmp.replace("a", ">").replace("b", "<");
            	writerRD.write(tmp + "\r\n");
            } else {
            	writerRD.write(textline + "\r\n");
            }

        }
        buffer.close();
        writerRD.close();
	}

}
