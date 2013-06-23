package uk.ac.ox.osscb.util;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;

public class TextUtil {
	
	public static void writeTextFile(String fileName, String txt){
		
	    Writer out = null; 
	    try {
	    	out = new OutputStreamWriter(new FileOutputStream(fileName));
			out.write(txt);
		} catch (IOException e) {
			throw new RuntimeException(String.format("Problem writing text of length: %d to a file: %s",
					txt.length(), fileName));
		} finally {
		      try {
				out.close();
			} catch (IOException e) {
				throw new RuntimeException(String.format(
						"Problem closing output stream writer when writing text of length: %d to a file: %s",
						txt.length(), fileName));
			}
	    }
	}

	/**
	 * TODO: close Scanner as well.
	 * @param path can be either relative or absolute
	 * @return
	 */
	public static String load(Class<?> clazz, String path){

		InputStream is = clazz.getResourceAsStream(path);

		if(null == is){
			String errMsg = String.format(
					"Could not find resource by path: '%s', loaded by class: %s",
					path, clazz.getName());
			throw new IllegalArgumentException(errMsg);
		}
		
		
		 try {
		        return new java.util.Scanner(is).useDelimiter("\\A").next();
		    } catch (java.util.NoSuchElementException e) {
		        return "";
		    }

		// Scanner sc = null;
//		try{
//			sc = new Scanner(is);
//			// String nextLine = sc.nextLine();
//			sc.useDelimiter("\\A");// read till the end
//			if(sc.hasNext()){
//				return sc.next();
//			}
//			
//		}finally{
//			try {
//				is.close();
//			} catch (IOException ex) {
//				throw new RuntimeException(String.format("Error closing stream. Path: '%s', class: %s.%n%s",
//						path, clazz.getName(), ex.getMessage()), ex);
//			}
//		}
//		
//		return path;
	}
}
