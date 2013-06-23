package uk.ac.ox.osscb.vienna;

import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.Util;

public class ViennaDotPsFileParser {
	
	private static final Logger log = LoggerFactory.getLogger(ViennaDotPsFileParser.class);
	
	private static final String beginLine = "%start of base pair probability data";
	
	// 0.32 1.00 hsb 1 39 0.042909 ubox
	private static final Pattern dataPattern = Pattern.compile("([\\d.]+)\\s+([\\d.]+)\\s+hsb\\s+(\\d+)\\s+(\\d+)\\s+([\\d.]+)\\s+(\\w+)");

	private Integer seqLength;
	
	/**
	 * Whether to subtract 1 from alidot.ps indices. Useful if we assume that
	 * those indices are one-based rather than zero-based.
	 */
	private final boolean decrementIndices = true;
	

	/**
	 * 
	 * @param seqLength can be null (dimensions of the output mtx are determined on the basis of vienna output. Not safe at all) 
	 */
	public ViennaDotPsFileParser(Integer seqLength) {
		super();
		this.seqLength = seqLength;
	}

	public double[][] parse(String viennaOutput){
		if(log.isDebugEnabled()){
			log.debug(String.format("Parsing vienna output:%n%s", viennaOutput));
		}
		Scanner scanner = new Scanner(viennaOutput);
		
		// use unix-style delimiter, windows will be matched
		scanner.useDelimiter("\n");
		// boolean isDataArea = false;

		List<ProbabilityPoint> probPoints = new LinkedList<ProbabilityPoint>();
		
		int[] maxIndices = new int[]{-1, -1};// 1st and 2nd resulting mtx correspondingly
		
		while(scanner.hasNext()){
			
			String line = scanner.next();
			
			ProbabilityPoint probabilityPoint = parseLine(line, maxIndices);
			log.debug("max indices on exit: {}, {}", maxIndices[0], maxIndices[1]);
			
			if(null != probabilityPoint)
				probPoints.add(probabilityPoint);

//			if(!isDataArea){
//				// we're in non-data area, check whether we got it
//				if(isDataArea = line.startsWith(beginLine)){
//					// we came accross data area, proceed to read the following lines
//					isDataArea = true;
//					continue;
//				}
//			}else{
//
//				ProbabilityPoint probabilityPoint = parseLine(line, maxIndices);
//				if(null != probabilityPoint)
//					probPoints.add(probabilityPoint);
//			}
		}
		log.debug("finished parsing. max Idx: " + maxIndices[0] + ", " +maxIndices[1]);
		
		int dim0, dim1;
		if(null != this.seqLength){
			dim0 = dim1 = this.seqLength;
			log.debug("using provided seq len: {}", this.seqLength);
		}else{
			dim0 = maxIndices[0]+1;
			dim1 = maxIndices[1]+1;
			log.debug("using parsed dimenstions: [{} x {}]", dim0, dim1);
		}
		
		double[][] probMtx = new double[dim0][dim1];
		
		// hack to avoid having zero probabilities
//		for(int i = 0; i < dim0; i++)
//			for(int j = 0; j < dim1; j++)
//				probMtx[i][j] = 1e-10;
		
		for(ProbabilityPoint dataPt: probPoints){
			// fill in the matrix
			probMtx[dataPt.getIndx0()][dataPt.getIndx1()] = dataPt.getVal();
		}

		if(log.isDebugEnabled()){
			log.debug(String.format("Got the following vienna probs mtx:%n%s", Util.print2DArray(probMtx)));
		}
		
		return probMtx;
	}

	private ProbabilityPoint parseLine(String line, int[] maxIndices) {
		
		if(maxIndices.length != 2){
			throw new IllegalArgumentException(String.format("Expected to get array of len 2, got: %d", maxIndices.length));
		}
		
		ProbabilityPoint probabilityPoint = null;
		
		Matcher dataM = dataPattern.matcher(line);
		
		if(dataM.matches()){
			log.debug("Line matches: {}", line);
			String group3 = dataM.group(3);
			String group4 = dataM.group(4);
			String group5 = dataM.group(5);
			String group6 = dataM.group(6).toLowerCase();
						
			if(group6.equals("ubox")){
				int idx0 = parseInt(group3, line);
				int idx1 = parseInt(group4, line);
				if(this.decrementIndices){
					idx0 = idx0 - 1;
					idx1 = idx1 - 1;
				}
				double dblVal = parseDbl(group5, line);
				
				if(idx0 > maxIndices[0])
					maxIndices[0] = idx0;
				if(idx1 > maxIndices[1])
					maxIndices[1] = idx1; 
				
				probabilityPoint = new ProbabilityPoint(idx0, idx1, dblVal);
				
			}else if(!group6.equals("lbox")){
				log.error(String.format("got a strange box type, omitting it: %s. Line: %s", group6, line));
			}			
		}else{
			log.debug("Line does NOT matches: {}", line);
		}

		log.debug("max indices {}, {}", maxIndices[0], maxIndices[1]);

		return probabilityPoint;
	}

	/**
	 * Probability with its coordinates
	 * 
	 * @author Vladimir
	 */
	private class ProbabilityPoint{
		private int indx0;
		private int indx1;
		private double val;
		public ProbabilityPoint(int indx0, int indx1, double val) {
			super();
			this.indx0 = indx0;
			this.indx1 = indx1;
			this.val = val;
		}
		public int getIndx0() {
			return indx0;
		}
		public int getIndx1() {
			return indx1;
		}
		public double getVal() {
			return val;
		}
	}

	public int parseInt(String raw, String context){
		try{
			int parsedInt = Integer.parseInt(raw);
			
			return parsedInt;
		}catch(NumberFormatException nfe){
			throw new IllegalArgumentException(String.format("Couldn't parse '%s' to int. Error: %s. Analysed line:%n\t%s",
					raw, nfe.getMessage(), context), nfe);
		}
	}

	public double parseDbl(String raw, String context){
		try{
			double parsedDbl = Double.parseDouble(raw);
			
			return parsedDbl;
		}catch(NumberFormatException nfe){
			throw new IllegalArgumentException(String.format(
					"Couldn't parse '%s' to double. Error: %s. Analysed line:%n\t%s",
					raw, nfe.getMessage(), context), nfe);
		}
	}
}
