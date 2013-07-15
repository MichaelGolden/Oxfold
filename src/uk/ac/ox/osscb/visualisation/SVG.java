package uk.ac.ox.osscb.visualisation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;

import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.image.JPEGTranscoder;
import org.apache.batik.transcoder.image.PNGTranscoder;

public class SVG {

	String bodyElement = "";

	double x = 0;
	double y = 0;
	double width = 0;
	double height = 0;

	String style = "";
	
	public SVG()
	{
		
	}
	
	public SVG clone()
	{
		SVG clone = new SVG();
		clone.x = this.x;
		clone.y = this.y;
		clone.width = this.width;
		clone.height = this.height;
		clone.bodyElement = this.bodyElement;
		clone.style = this.style;
		
		return clone;
	}

	public SVG(double width, double height) {
		this.width = width;
		this.height = height;
	}

	public String getSimpleStartElement() {
		return "<svg x=\"" + x + "\" y=\"" + y + "\" style=\"" + style + "\">";
	}

	private String getStartElement() {
		return "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" x=\""
				+ x + "\" y=\"" + y + "\" width=\"" + width + "\" height=\""
				+ height + "\"  style=\"" + style + "\">";
	}

	public void setMainStyle(String style) {
		this.style = style;
	}

	public void setBodyElement(String bodyElement) {
		this.bodyElement = bodyElement;
	}

	public void setDimensions(double width, double height) {
		this.width = width;
		this.height = height;
	}

	public void setLocation(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public String getSimpleSVGString() {
		return getSimpleStartElement() + "\n" + bodyElement + "</svg>";
	}

	public String getSVGString() {
		return getStartElement() + "\n" + bodyElement + "</svg>";
	}

	public void writeSVG(File svgFile) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(svgFile));
		writer.write(getSVGString());
		writer.close();
	}
	
	public SVG attachRight(SVG attachment, double spacer)
	{
		SVG ret = new SVG(this.width+spacer+attachment.width, Math.max(this.height, attachment.height));

		SVG clone = attachment.clone();
		clone.x = this.width+spacer;
		
		ret.bodyElement = this.getSimpleSVGString()+"\n\n"+clone.getSimpleSVGString();

		return ret;
	}
	
	public SVG attachBottom(SVG attachment, double spacer)
	{
		SVG ret = new SVG(Math.max(this.width, attachment.width), this.height+spacer+attachment.height);
		//SVG ret = new SVG(this.width+spacer+attachment.width, Math.max(this.height, attachment.height));

		SVG clone = attachment.clone();
		clone.y = this.height+spacer;
		
		ret.bodyElement = this.getSimpleSVGString()+"\n\n"+clone.getSimpleSVGString();

		return ret;
	}


	/*
	public void savePNG(File svgFile, File pngFile) throws IOException {
		writeSVG(svgFile);

		
	
		//String[] cmds = {"binaries/", "/c",""+INKSCAPE_EXECUTABLE_PATH+"", "-z", "-e", "\""+svgFile.getAbsolutePath()+"\"", "\""+pngFile.getAbsolutePath()+"\""};
		
		//Process process = new ProcessBuilder(cmds).start();
		//System.
		//ProcessBuilder pb = new ProcessBuilder().
		String cmd = "java -jar binaries/batik-1.7/batik-rasterizer.jar \"" +svgFile.getAbsolutePath() + "\"";
		System.out.println(cmd);
		Process process = Runtime.getRuntime().exec(cmd);
		BufferedReader r = new BufferedReader(new InputStreamReader(process.getErrorStream()));
        String textline = null;
        while((textline = r.readLine()) != null)
        {
        	System.out.println("error: "+textline);
        }
		try {
			int exitCode = process.waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}*/
	
	public void savePNG(File svgFile, File pngFile) throws IOException {
		writeSVG(svgFile);
		
	
		PNGTranscoder t = new PNGTranscoder();
		//t.addTranscodingHint(PNGTranscoder.KEY_PIXEL_UNIT_TO_MILLIMETER, new Float(0.1));	
		t.addTranscodingHint(PNGTranscoder.KEY_WIDTH, new Float(Math.min(this.width, 10000)));

        // Create the transcoder input.
        String svgURI = svgFile.toURL().toString();
        TranscoderInput input = new TranscoderInput(svgURI);

        // Create the transcoder output.
        OutputStream ostream = new FileOutputStream(pngFile);
        TranscoderOutput output = new TranscoderOutput(ostream);

        // Save the image.
        try {
			t.transcode(input, output);
		} catch (TranscoderException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

        // Flush and close the stream.
        ostream.flush();
        ostream.close();
	}
}
