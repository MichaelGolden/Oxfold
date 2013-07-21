package uk.ac.ox.osscb.visualisation;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;

import uk.ac.ox.osscb.analysis.BasePairMetrics;
import uk.ac.ox.osscb.analysis.MountainMetrics;
import uk.ac.ox.osscb.analysis.RNAFoldingTools;
import uk.ac.ox.osscb.analysis.StructureData;

public class DataVisualiser {
	
	public static final ColorGradient alphaRedGradient =  new ColorGradient(new Color(255,255,255,0), new Color(255,0,0,255));
	public static final ColorGradient whiteRedGradient =  new ColorGradient(new Color(255,255,255,255), new Color(255,0,0,255));
	
	Color missingDataColor = Color.gray;
    Color filteredDataColor = Color.darkGray;
	Color nucleotideBackgroundColor = missingDataColor;
	
    Color bondColor = Color.gray;
    int bondThickness = 4;
    
    int nucleotideDiameter = 40;
    
    boolean showDNA = false;
    NucleotideComposition.Type nucleotideCompositionType =  NucleotideComposition.Type.FREQUENCY;	
    
    double xoffset = 150;
   // double zoomScale = 1;
    double horizontalScale = 2.6;
    double verticalScale = 2.6;
    
	
	double minx;
	double miny;
	double maxx;
	double maxy;
	
	public enum DrawingMode {

        NAVIEW, RADIATE_VIEW, RADIATE_VIEW_FLAT
    };
    DrawingMode drawingMode = DrawingMode.RADIATE_VIEW;
	
	public Point2D.Double [] computeStructureToBeDrawn(int[] pairedSites) {
		ArrayList<Point2D.Double> np = null;
		Point2D.Double [] nucleotidePositions = null;

		
        try {
            switch (drawingMode) {
                case NAVIEW:
                    np = NAView.naview_xy_coordinates(pairedSites);
                    break;
                case RADIATE_VIEW:
                    np = RadiateView.radiateview_xy_coordinates(pairedSites, false);
                    break;
                case RADIATE_VIEW_FLAT:
                    np = RadiateView.radiateview_xy_coordinates(pairedSites, true);
                    break;
            }
            // 
        } catch (Exception ex) {
            Logger.getLogger(SubstructureDrawPanel.class.getName()).log(Level.SEVERE, null, ex);
        }

        minx = Double.MAX_VALUE;
        miny = Double.MAX_VALUE;
        maxx = Double.MIN_VALUE;
        maxy = Double.MIN_VALUE;

        for (int i = 0; i < np.size(); i++) {
            Point2D.Double pos = np.get(i);
            minx = Math.min(minx, pos.x);
            miny = Math.min(miny, pos.y);
            maxx = Math.max(maxx, pos.x);
            maxy = Math.max(maxy, pos.y);
        }
        nucleotidePositions = new Point2D.Double[np.size()];
        for (int i = 0; i < nucleotidePositions.length; i++) {
            nucleotidePositions[i] = new Point2D.Double();
            nucleotidePositions[i].x = xoffset + (np.get(i).x - minx) * horizontalScale;
            nucleotidePositions[i].y = 50 + (np.get(i).y - miny) * verticalScale;
        }       
        
        return nucleotidePositions;
    }
	
	public static SVG drawBasePairProbMatrix(double [][] basePairProb, double width, double height)
	{
		SVG matrix = new SVG(width, height);
		
		double cellWidth = width / ((double)basePairProb.length);
		double cellHeight = height / ((double)basePairProb[0].length);
		
		StringBuffer sb = new StringBuffer();
		sb.append("<rect x=\"0\" y=\"0\" width=\""+width+"\" height=\""+height+"\" style=\"fill:rgb(200,200,200)\"/>");
		sb.append("\n");
		matrix.bodyElement += "";
		
		for(int i = 0 ; i < basePairProb.length ; i++)
		{
			for(int j = 0 ; j < basePairProb[0].length ; j++)
			{
				double x = i*cellWidth;
				double y = j*cellHeight;
				Color c = alphaRedGradient.getColor((float)basePairProb[i][j]);
				sb.append("<rect x=\"" + x+"\" y=\""+y+"\" width=\""+cellWidth+"\" height=\""+cellHeight+"\" style=\"fill:" + GraphicsUtils.getRGBAString(c) + "\"/>\n");				
			}
			//System.out.println(i+"\t"+basePairProb[0].length);
		}
		matrix.bodyElement = sb.toString();
		return matrix;
		
	}
	
	public SVG drawStructure(StructureData data, double [][] basePairProb)
	{
		int [] unpairedSites = new int[data.pairedSites.length];
		
		Point2D.Double [] nucleotidePositions = computeStructureToBeDrawn(unpairedSites);
		//Point2D.Double [] nucleotidePositions = computeStructureToBeDrawn(data.pairedSites);
		
		NucleotideComposition nucleotideComposition = null;
		if(data.sequences != null)
		{
			nucleotideComposition = new NucleotideComposition(data.sequences, data.sequenceNames);
		}
		double [] nucleotidePairedProbs = null;
		if(basePairProb != null)
		{
			nucleotidePairedProbs = StructureData.getNucleotidePairedProbs(basePairProb);
		}
		return drawStructure(nucleotidePositions, data.pairedSites, nucleotideComposition, null, true, true, basePairProb, nucleotidePairedProbs, minx, miny, maxx, maxy);
	}
	
	public SVG drawStructure(Point2D.Double [] nucleotidePositions, int [] pairedSites, NucleotideComposition nucleotideSource, String title, boolean drawSequenceLogo, boolean showBonds, double [][] basePairProb, double [] nucleotidePairingProbs, double minx, double miny, double maxx, double maxy) {

		
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);

        /*
        if (substructureModel.substructure == null || substructureModel.structureOverlay == null || substructureModel.structureOverlay.pairedSites == null || nucleotidePositions == null) {
            return null;
        }*/

        int panelWidth = (int) ((maxx - minx) * horizontalScale + xoffset * 2);
        int panelHeight = (int) ((maxy - miny) * verticalScale + 100);


		SVG svg = new SVG(panelWidth, panelHeight);
		svg.setMainStyle("fill:none;stroke-width:16");

		if(title != null)
		{
			pw.println("<g>");
			pw.print("<text font-family=\"sans-serif\" font-size=\"16\" x=\""+(panelWidth/2)+"\" y=\""+16+"\" style=\"text-anchor: middle; fill:rgb(0,0,0);\">");
			pw.print(title);
			pw.println("</text>");
			pw.println("</g>");
		}
        
        // draw two-dimensional data
        pw.println("<g>");
        int length = pairedSites.length;
        if (basePairProb != null) {
            for (int k = 0 ; k < pairedSites.length ; k++) {
                for (int l = k ; l < pairedSites.length ; l++) {

                   // if (substructureModel.maxDistance == -1 || (substructureModel.substructureDistanceMatrix != null && substructureModel.substructureDistanceMatrix.getDistance(k, l) <= substructureModel.maxDistance) || (substructureModel.substructureDistanceMatrix == null && substructureModel.substructureDistanceMatrix.getDistance(i, j) <= substructureModel.maxDistance)) {
                        Color c = null;

                        double value = basePairProb[k][l];
                       // double value = substructureModel.data2D.get(i, j, substructureModel.mapping2D);
                        if (value < 1e-4) {
                            c = null;
                        } else if (value > 0) {
                           c = alphaRedGradient.getColor((float)value);
                        }

                        if (c != null) { 

                            double x1 = nucleotidePositions[k].getX();
                            double y1 = nucleotidePositions[k].getY();
                            double x2 = nucleotidePositions[l].getX();
                            double y2 = nucleotidePositions[l].getY();

                            int structureMidpoint = pairedSites.length / 2;
                            pw.println("    <polyline id=\"2d_data_" + k + "_" + l + "\" points=\"" + x1 + " " + y1 + " " + x2 + " " + y2 + "\" style=\"stroke-width:7.5;stroke:" + GraphicsUtils.getRGBAString(c) + "\"/>");
                        }
                    }
            }
        }
        pw.println("</g>");
        
		
        if (showBonds) {
            pw.println("<g>");
            for (int i = 0; i < nucleotidePositions.length; i++) {
            	  int a = i;
                  int b = pairedSites[i] - 1;
                  if (i + 1 < pairedSites[i]) {
                      /*
                       * Line2D bond = new Line2D.Double(nucleotidePositions[a],
                       * nucleotidePositions[b]); g.setStroke(new BasicStroke(2));
                       * g.setColor(Color.gray); g.draw(bond);
                       */
                      pw.println("    <line id=\"bond_" + i + "\" x1=\"" + nucleotidePositions[a].x + "\" y1=\"" + nucleotidePositions[a].y + "\"  x2=\"" + nucleotidePositions[b].x + "\" y2=\"" + nucleotidePositions[b].y + "\" style=\"stroke-width:" + bondThickness + ";stroke:" + GraphicsUtils.getRGBAString(bondColor) + "\"/>");
                      //pw.println("    <polyline id=\"bond_" + i + "\" points=\"" + nucleotidePositions[a].x + " " + nucleotidePositions[a].y + " " + nucleotidePositions[b].x + " " + nucleotidePositions[b].y + "\" style=\"stroke-width:\"" + bondThickness + "\";stroke:#" + getHexString(bondColor) + "\"/>");
                  }
            }
            pw.println("</g>");
        }


        // draw the nucleotides
        pw.println("<g>");
        for (int i = 0; i < nucleotidePositions.length; i++) {
            Color nucleotideBackgroundColor = missingDataColor;
        	if(nucleotidePairingProbs != null)
        	{
        		double p = nucleotidePairingProbs[i];
        		nucleotideBackgroundColor = whiteRedGradient.getColor((float)p);       

        	}
        	

    		pw.println("    <circle id=\"nucleotide_" + (i) + "\" cx=\"" + nucleotidePositions[i].getX() + "\" cy=\"" + nucleotidePositions[i].getY() + "\" r=\"" + nucleotideDiameter / 2 + "\" style=\"stroke-width:2;stroke:black;fill:" + GraphicsUtils.getRGBAString(nucleotideBackgroundColor) + "\"/>");

            // draw the information
            /*
             * if (drawSequenceLogo) { Color bestColor =
             * ColorUtils.selectBestForegroundColor(nucleotideBackgroundColor,
             * Color.white, Color.black); if (model.nucleotideComposition !=
             * null) { if (model.nucleotideCompositionType ==
             * NucleotideComposition.Type.SHANNON_ENTROPY) { // fa =
             * structure.shannonFrequencies[i]; double[] fa =
             * Arrays.copyOf(model.nucleotideComposition.mappedShannonComposition[(model.structure.startPosition
             * + i - 1) % model.sequenceLength], 5); for (int k = 0; k < 4; k++)
             * { fa[k] = fa[k] / 2; } pw.print(drawSequenceLogoSVG("logo_" +
             * (model.structure.startPosition + i),
             * nucleotidePositions[i].getX(), nucleotidePositions[i].getY() -
             * (nucleotideDiameter / 2) + 4, nucleotideDiameter,
             * nucleotideDiameter - 8, fa, bestColor)); } else if
             * (model.nucleotideCompositionType ==
             * NucleotideComposition.Type.FREQUENCY) { double[] fa =
             * model.nucleotideComposition.mappedFrequencyComposition[(model.structure.startPosition
             * + i - 1) % model.sequenceLength];
             * pw.print(drawSequenceLogoSVG("logo_" +
             * (model.structure.startPosition + i),
             * nucleotidePositions[i].getX(), nucleotidePositions[i].getY() -
             * (nucleotideDiameter / 2) + 4, nucleotideDiameter,
             * nucleotideDiameter - 8, fa, bestColor)); } } }
             */

            
            if (drawSequenceLogo) {
                Color bestColor = ColorUtils.selectBestForegroundColor(nucleotideBackgroundColor, Color.white, Color.black);
                int nucPos = i;
                if (nucleotideSource != null) {
                    if (nucleotideCompositionType == NucleotideComposition.Type.SHANNON_ENTROPY) {
                        // fa = structure.shannonFrequencies[i];
                        double[] nucfa = nucleotideSource.shannonComposition[i];
                        if (nucfa != null) {
                            double[] fa = Arrays.copyOf(nucfa, 5);
                            for (int k = 0; k < 4; k++) {
                                fa[k] = fa[k] / 2;
                            }
                            pw.print(drawSequenceLogoSVG("logo_" + (i), nucleotidePositions[i].getX(), nucleotidePositions[i].getY() - (nucleotideDiameter / 2) + 3, nucleotideDiameter, nucleotideDiameter - 8, fa, bestColor));
                        }
                    } else if (nucleotideCompositionType == NucleotideComposition.Type.FREQUENCY) {
                        double[] nucfa =nucleotideSource.frequencyComposition[i];
                        if (nucfa != null) {
                            pw.print(drawSequenceLogoSVG("logo_" + (i), nucleotidePositions[i].getX(), nucleotidePositions[i].getY() - (nucleotideDiameter / 2) + 3, nucleotideDiameter, nucleotideDiameter - 8, nucfa, bestColor));
                        }
                    }
                }
            }
        }
        pw.println("</g>");

        /*
        // draw position lines
        pw.println("<g>");
        for (int i = 0; i < nucleotidePositions.length; i++) {
            int offsetx = 0;
            double side = 1;
            String textanchor = "start";
            if (i < length / 2) {
                offsetx = -(int) (nucleotideDiameter / 2) - 5;
                side = -1;
                textanchor = "end";
            } else {
                offsetx = (int) (nucleotideDiameter / 2) + 5;
            }

            if (nucleotidePositions[i] != null) {
                int pos = (substructureModel.substructure.getStartPosition() + i) % substructureModel.sequenceLength + 1;
               // int pos = (substructureModel.substructure.getStartPosition() + i - 1) % substructureModel.sequenceLength + 2;
                double fontSize = 11;
                if (substructureModel.numbering != 0 && pos % substructureModel.numbering == 0) {
                    pw.println("    <text id=\"nucleotide_position_" + pos + "\" x=\"" + (offsetx + nucleotidePositions[i].getX()) + "\" y=\"" + (nucleotidePositions[i].getY() + (fontSize / 2)) + "\" style=\"font-size:" + fontSize + "px;stroke:none;fill:black\" text-anchor=\"" + textanchor + "\" >");
                    pw.println("        <tspan>" + pos + "</tspan>");
                    pw.println("    </text>");

                    double x1 = nucleotidePositions[i].getX() + (side * nucleotideDiameter / 2) - 2.5;
                    double y1 = nucleotidePositions[i].getY();
                    double x2 = nucleotidePositions[i].getX() + (side * nucleotideDiameter / 2) + 2.5;
                    double y2 = nucleotidePositions[i].getY();
                    pw.println("    <polyline points=\"" + x1 + " " + y1 + " " + x2 + " " + y2 + "\" style=\"stroke-width:1;stroke:black\"/>");

                }
            }
        }
        pw.println("</g>");
        */
        
        pw.close();
        svg.bodyElement = sw.toString();
        return svg;
    }
	
	public String drawSequenceLogoSVG(String id, double x, double y, double width, double height, double[] h, Color textColor) {
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);

        double fontHeight = 32;
        double scale = (height / fontHeight);
        double base = y;
        for (int i = 0; i < h.length; i++) {
            //h[i] = 0.25;
            double fontHeightScale = (h[i]);
            //Font tallerFont = f1.deriveFont(AffineTransform.getScaleInstance(scale, fontHeightScale * scale));

            String a = "";
            switch (i) {
                case 0:
                    a = "A";
                    break;
                case 1:
                    a = "C";
                    break;
                case 2:
                    a = "G";
                    break;
                case 3:
                    if (showDNA) {
                        a = "T";
                    } else {
                        a = "U";
                    }
                    break;
            }
            String b = a.length() > 0 ? a : "X";
            
            base += fontHeightScale * scale * fontHeight;

            if (!b.equals("X") && h[i] > 0) {
                float xf = (float) x;
                float yf = (float) (base);
                xf /= scale;
                yf /= fontHeightScale * scale;
                //  pw.println("    <g transform=\"scale(" + scale + "," + fontHeightScale * scale + ")\">");
                pw.print("        <text transform=\"scale(" + scale + "," + fontHeightScale * scale * 1.1 + ")\" x=\"" + xf + "\" y=\"" + (yf / 1.1) + "\" id=\"" + id + "_" + i + "\"  style=\"font-size:" + (int) (fontHeight) + ";stroke:none;fill:" + GraphicsUtils.getRGBAString(textColor) + ";text-anchor:middle;\">");
                pw.println("            <tspan id=\"base\">" + b + "</tspan>");
                pw.println("</text>");
            }
        }

        return sw.toString();
    }
	
	public static String getMetricString(StructureData predicted, StructureData experimental)
	{

		double sensitivity = BasePairMetrics.calculateSensitivity(experimental.pairedSites, predicted.pairedSites);
		double ppv = BasePairMetrics.calculatePPV(experimental.pairedSites, predicted.pairedSites);
		double fscore = BasePairMetrics.calculateFScore(experimental.pairedSites,predicted.pairedSites);
		double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(experimental.pairedSites, predicted.pairedSites);
		String 	metrics = "Sensitivity = "+sensitivity;
				metrics+= "\nPPV = "+ppv;
				metrics+= "\nFScore = "+fscore;
				metrics+= "\nMountain = "+mountainSim;
	   return metrics;
	}
	
	public SVG drawComparisonPredictedExperimental(StructureData predicted, StructureData experimental)
	{		
		double[][] basePairProb = null;
		
		
		try {
			basePairProb = predicted.basePairProbFile == null ? null : StructureData.getBasePairProb(predicted.basePairProbFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		SVG structure1 = drawStructure(predicted, basePairProb);
		SVG title1 = drawText(predicted.title, 24f, structure1.width/2, 24, Align.CENTER, "normal", false);
		structure1 = title1.attachBottom(structure1, 0);
		
		SVG structure2 = drawStructure(experimental, basePairProb);
		SVG title2 = drawText(experimental.title, 24f, structure2.width/2, 24, Align.CENTER,  "normal", false);
		structure2 = title2.attachBottom(structure2, 0);
		
		if(basePairProb != null)
		{
			//SVG basePairSVG = this.drawBasePairProbMatrix(basePairProb, structure2.width,  structure2.width);
			//structure2 = structure2.attachBottom(basePairSVG, 0);
		}
		

		SVG full = structure1.attachRight(structure2, 50);

		
		full = full.attachBottom(drawText(getMetricString(predicted, experimental), 24f, 20, 0, Align.LEFT,  "normal", false), 10);
		
		String structures = RNAFoldingTools.getDotBracketStringFromPairedSites(predicted.pairedSites)+"\n";
		structures += RNAFoldingTools.getDotBracketStringFromPairedSites(experimental.pairedSites)+"\n";		
		full = full.attachBottom(drawText(structures, 24f, 20, 0, Align.LEFT,  "0.3em", true), 10);
		return full;
	}
	
	public SVG drawComparisonPredicted(StructureData d1, StructureData d2, StructureData experimental)
	{		
		double[][] basePairProb1 = null;
		double[][] basePairProb2 = null;
		
		try {
			basePairProb1 = d1.basePairProbFile == null ? null : StructureData.getBasePairProb(d1.basePairProbFile);
			basePairProb2 = d2.basePairProbFile == null ? null : StructureData.getBasePairProb(d2.basePairProbFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		SVG structure1 = drawStructure(d1, basePairProb1);
		SVG title1 = drawText(d1.title, 24f, structure1.width/2, 24, Align.CENTER, "normal", false);
		structure1 = title1.attachBottom(structure1, 0);
		structure1 = structure1.attachBottom(drawText(getMetricString(d1, experimental), 24f, 20, 0, Align.LEFT,  "normal", false), 10);
		
		
		SVG structure2 = drawStructure(d2, basePairProb2);
		SVG title2 = drawText(d2.title, 24f, structure2.width/2, 24, Align.CENTER,  "normal", false);
		structure2 = title2.attachBottom(structure2, 0);
		structure2 = structure2.attachBottom(drawText(getMetricString(d2, experimental), 24f, 20, 0, Align.LEFT,  "normal", false), 10);
		
		SVG structure3 = drawStructure(experimental, null);
		SVG title3 = drawText(experimental.title, 24f, structure2.width/2, 24, Align.CENTER,  "normal", false);
		structure3 = title3.attachBottom(structure3, 0);
		//structure3 = structure3.attachBottom("", 24f, 20, 0, Align.LEFT,  "normal", false), 10);
		
		if(basePairProb1 != null && basePairProb2 != null);
		{
			//SVG basePairSVG = this.drawBasePairProbMatrix(basePairProb, structure2.width,  structure2.width);
			//structure2 = structure2.attachBottom(basePairSVG, 0);
		}
		

		SVG full = structure1.attachRight(structure2, 50).attachRight(structure3, 50);

		
		//full = full.attachBottom(drawText(text, 24f, 20, 0, Align.LEFT,  "normal", false), 10);
		
		String structures = RNAFoldingTools.getDotBracketStringFromPairedSites(d1.pairedSites)+"\n";
		structures += RNAFoldingTools.getDotBracketStringFromPairedSites(d2.pairedSites)+"\n";
		structures += RNAFoldingTools.getDotBracketStringFromPairedSites(experimental.pairedSites)+"\n";		
		full = full.attachBottom(drawText(structures, 24f, 20, 0, Align.LEFT,  "0.3em", true), 10);
		return full;
	}
	
	enum Align{
		CENTER("middle"), LEFT("start"), RIGHT("end");
		
		String textAnchor;
		
		Align(String textAnchor)
		{
			this.textAnchor = textAnchor;
		}
		
		public String toString()
		{
			return this.textAnchor;
		}
	}
	
	public SVG drawText(String text, float fontSize, double x, double y, Align textAlignment, String letterSpacing, boolean monospaceFont)
	{		
		String [] lines = text.split("(\n)+");
		
		//double height = fontSize*lines.length;
		
		SVG textSVG = new SVG();
		double lineSpacing = fontSize*0.2;
		textSVG.height = (fontSize+lineSpacing)*(lines.length+1);
		
		String fontFamily = "sans-serif";
		if(monospaceFont)
		{
			fontFamily = "monospace";
		}
		
		textSVG.bodyElement += "<text font-family=\""+fontFamily+"\" font-size=\""+fontSize+"\" x=\""+x+"\" y=\""+y+"\" style=\"text-anchor: "+textAlignment.toString()+"; fill:rgb(0,0,0) ; letter-spacing:"+letterSpacing+";\">\n";
		for(int i = 0 ; i < lines.length ; i++)
		{
			//textSVG.bodyElement += "<tspan x=\""+x+"\" y=\""+fontSize*i+"\">"+lines[i]+"</tspan>\n";
			textSVG.bodyElement += "<tspan x=\""+x+"\" dy=\""+(fontSize+lineSpacing)+"\">"+lines[i]+"</tspan>\n";
		}
		textSVG.bodyElement += "</text>\n";
		
		return textSVG;
		
	}
	
	
	
	public static void main(String [] args)
	{
		String structure = "(....)(((((((((((((((((.........)((((((..(((.............)).....))))))).............))))))(.)(((((((((((..............)))))))))))..)((((((((.......))))))))))))(((((.............................))))).....)))))......";
		DataVisualiser visualiser = new DataVisualiser();
		//visualiser.computeStructureToBeDrawn(RNAFoldingTools.getPairedSitesFromDotBracketString(structure));
		
		//BufferedWriter writer = new BufferedWriter(new FileWriter(new File("test.svg")));
		//writer.wr
		
		File dataDir = new File("datasets/");
		File outputDir = new File("output3/");
		File resultsFile = new File("resultsCofoldWeighted.csv");
		outputDir.mkdir();
		System.out.println(dataDir.list().length);
		ArrayList<StructureData> experimentalStructures = new ArrayList<StructureData>();
		for(File experimentalFile : dataDir.listFiles())
		{
			try {
				experimentalStructures.add(StructureData.readExperimentalStructureData(experimentalFile));				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		Collections.sort(experimentalStructures);
		
		//SVG full = visualiser.drawComparison(experimentalStructures.get(0), experimentalStructures.get(1));
		/*try {
			full.savePNG(new File("test.svg"), new File("test.png"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		//System.out.println(full.getSVGString());
		/*
		for(int i = 0 ; i < experimentalStructures.size() ; i++)
		{
			StructureData data = experimentalStructures.get(i);
			System.out.println(visualiser.drawStructure(data));
			break;
		}		*/
	}
}
