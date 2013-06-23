package uk.ac.ox.osscb;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ox.osscb.util.ProbabilityValueValidator;

public class InsideOutsideProbabilities {

	private final Logger log = LoggerFactory.getLogger(InsideOutsideProbabilities.class);

	private Map<Character, BigDecimal[][]> nonTerminalIndices;
	private final int dimension;
	
	private MathContext mathCtx = MathContext.DECIMAL32; //null;// new MathContext(5);
	
	public String testToStr(){
		BigDecimal[][] bigDecimals = this.nonTerminalIndices.get('S');
		return bigDecimals[0][bigDecimals[0].length-1].toString();
	}
	
	public String printAllTables(){
		StringBuilder sb = new StringBuilder();
		for(Map.Entry<Character, BigDecimal[][]> entry : nonTerminalIndices.entrySet()){
			sb.append(Util.nL()).append( entry.getKey()).append(":");
			BigDecimal[][] numAr = entry.getValue();
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
		
		this.nonTerminalIndices = new HashMap<Character, BigDecimal[][]>();
		for(char c : nonTerminals){
			if(this.nonTerminalIndices.containsKey(c)){
				throw new RuntimeException("TODO-3");
			}
			
			this.nonTerminalIndices.put(c, Util.makeSquareZerosMatrix(strLen));
		}		

		this.dimension = strLen;
	}

	public void setProb(char nonTerminal, int i, int j, double prob){
		setProb(nonTerminal, i, j, BigDecimal.valueOf(prob));
	}
	
	public void setProb(char nonTerminal, int i, int j, BigDecimal prob){
		ProbabilityValueValidator.validateP(prob);
		//log.debug("'%c'[%d, %d]->%.5f", new Object[]{nonTerminal, i, j, prob});
		log.debug("'{}'[{}, {}]->{}", new Object[]{nonTerminal, i, j, prob});
		assertNonTerminalExists(nonTerminal, i, j);
		
		BigDecimal[][] probs = this.nonTerminalIndices.get(nonTerminal);
		probs[i][j] = null != this.mathCtx ? prob.round(this.mathCtx) : prob;
	}
	
	public BigDecimal getProb(char nonTerminal, int i, int j){
		assertNonTerminalExists(nonTerminal, i,j );
		
		BigDecimal p = this.nonTerminalIndices.get(nonTerminal)[i][j];
		// ProbabilityValueValidator.validateP(p);
		return p;
	}
	
	public void increment(char nonTerminal, int i, int j, BigDecimal incrementVal){
		assertNonTerminalExists(nonTerminal, i,j );
		
		BigDecimal initialProb = getProb(nonTerminal, i, j);
		BigDecimal newProb = initialProb.add(incrementVal);
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
}
