package uk.ac.ox.osscb;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.util.ProbabilityValueValidator;

public class InsideOutsideProbabilities {

	private final Logger log = LoggerFactory.getLogger(InsideOutsideProbabilities.class);

	private Map<Character, PointRes[][]> nonTerminalIndices;
	private final int dimension;
	
	private MathContext mathCtx = MathContext.DECIMAL32; //null;// new MathContext(5);
	
	public String testToStr(){
		PointRes[][] bigDecimals = this.nonTerminalIndices.get('S');
		return bigDecimals[0][bigDecimals[0].length-1].toString();
	}
	
	public String printAllTables(){
		StringBuilder sb = new StringBuilder();
		for(Map.Entry<Character, PointRes[][]> entry : nonTerminalIndices.entrySet()){
			sb.append(Util.nL()).append( entry.getKey()).append(":");
			PointRes[][] numAr = entry.getValue();
			sb.append(Util.print2DArray(numAr));
		}
		return sb.toString();
	}

	public InsideOutsideProbabilities(char[] nonTerminals, int strLen){
		if(null == nonTerminals || nonTerminals.length < 1){
			throw new RuntimeException("TODO-1");
		}
		
		if(strLen < 0){
			throw new RuntimeException("TODO-2");
		}
		
		this.nonTerminalIndices = new HashMap<Character, PointRes[][]>();
		for(char c : nonTerminals){
			if(this.nonTerminalIndices.containsKey(c)){
				throw new RuntimeException("TODO-3");
			}
			
			this.nonTerminalIndices.put(c, Util.makeSquareZerosMatrix(strLen));
		}		

		this.dimension = strLen;
	}

	public void setProb(char nonTerminal, int i, int j, double prob){
		setProb(nonTerminal, i, j, PointRes.valueOf(prob));
	}
	
	public void setProb(char nonTerminal, int i, int j, PointRes prob){
		ProbabilityValueValidator.validateP(prob);
		//log.debug("'%c'[%d, %d]->%.5f", new Object[]{nonTerminal, i, j, prob});
		log.debug("'{}'[{}, {}]->{}", new Object[]{nonTerminal, i, j, prob});
		assertNonTerminalExists(nonTerminal, i, j);
		
		PointRes[][] probs = this.nonTerminalIndices.get(nonTerminal);
		probs[i][j] = null != this.mathCtx ? prob.round(this.mathCtx) : prob;
	}
	
	public PointRes getProb(char nonTerminal, int i, int j){
		assertNonTerminalExists(nonTerminal, i,j );
		
		PointRes p = this.nonTerminalIndices.get(nonTerminal)[i][j];
		// ProbabilityValueValidator.validateP(p);
		return p;
	}
	
	public void increment(char nonTerminal, int i, int j, PointRes incrementVal){
		assertNonTerminalExists(nonTerminal, i,j );
		
		PointRes initialProb = getProb(nonTerminal, i, j);
		PointRes newProb = initialProb.add(incrementVal);
		setProb(nonTerminal, i, j, newProb);
	}
	
	private void assertNonTerminalExists(char nonTerminal, int i, int j){

		if(!this.nonTerminalIndices.containsKey(nonTerminal)){
			throw new RuntimeException(String.format("non-terminal %s is not found. Available are: [%s]",
					nonTerminal, StringUtils.join(this.nonTerminalIndices.keySet(), ";\t")));  
		}

		if(i >= this.dimension){
			throw new RuntimeException(String.format("i: %d is bigger than dimenstions allowed: %d", i, this.dimension));
		}
		if(j >= this.dimension){
			throw new RuntimeException(String.format("j: %d is bigger than dimenstions allowed: %d", j, this.dimension));
		}
	}

	public int getDimension() {
		return dimension;
	}
	
	public void writeTable(File outFile, char nonTerminal) throws IOException
	{
		DecimalFormat df = new DecimalFormat("0.000E000");
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		PointRes [][] matrix = nonTerminalIndices.get(new Character(nonTerminal));
		for(int i = 0 ; i < matrix.length ; i++)
		{
			for(int j = 0 ; j < matrix[0].length ; j++)
			{
				writer.write(df.format(matrix[i][j].doubleValue())+"\t");
			}
			writer.newLine();
		}
		writer.close();
	}
}
