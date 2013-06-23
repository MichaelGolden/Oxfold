import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


public class FastTree {
	public static String FAST_TREE_EXECUTABLE = "binaries/FastTree.exe";
	
	public static void nucleotideGTRFastTree(File fastaFile, File outNewick) throws Exception
	{
		  String cmd = "cmd /c " + new File(FAST_TREE_EXECUTABLE).getAbsolutePath() + " -nosupport -nt -gtr " + fastaFile.getAbsolutePath() + " > " + outNewick.getAbsolutePath();
          Process process = Runtime.getRuntime().exec(cmd);
          nullOutput(process.getInputStream());
          nullOutput(process.getErrorStream());
          int exitCode = process.waitFor();
          
	}
	
	   public static void nullOutput(final InputStream inputStream) {
	        new Thread() {

	            @Override
	            public void run() {
	                try {

	                    BufferedReader buffer = new BufferedReader(new InputStreamReader(inputStream));
	                    String textline = null;
	                    while ((textline = buffer.readLine()) != null) {
	                    }
	                    buffer.close();
	                } catch (IOException ex) {
	                }
	            }
	        }.start();
	    }
}
