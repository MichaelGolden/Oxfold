package uk.ac.ox.osscb;



import uk.ac.ox.osscb.grammar.RuleType;

/**
 * Production rule class
 * Class to represent production rules for a stochastic context-free grammar to use for 
 * RNA structures. 
 * 
 * Each rule is one of 3 classes:
 * 1: A -> . 
 * 2: A -> (B)
 * 3: A -> BC
 * where terminals are {(,.,)} and nonterminals are {A, B, C}
 * 16 July 2012
 * @author lmath@cs.ubc.ca
 *
 */

public class ProductionRule {

	/* For a sample rule A-> BC
	 * left = ['A'], right = ['B', 'C']
	 */
	char left;
	char[] right;
	RuleType ruleclass;
	PointRes probability;
	
	
	/**
	 * Constructor
	 * Input left, right as specified above and the associated probability for the rule.
	 */
	public ProductionRule(char left, char[] right, PointRes probability)
	{
		this.left = left;
		this.right = right;
		this.probability = probability;
		ruleclass = null;
	}
	
	/**
	 * Getter for the class of the production rule as explained at head of this file.
	 * @return type
	 */
	public RuleType getRuleClass()
	{
		return ruleclass;
	}
	
	public PointRes getProbability()
	{
		return probability;
	}
	
	public char getLeft()
	{
		return left;
	}
	
	public char[] getRight()
	{
		return right;
	}
	
	public void setRuleClass(RuleType ruleclass)
	{
		this.ruleclass = ruleclass;
	}
	

	/* Just wanted to see if Arrays.toString would print char arrays nicely
	 * and might want this syntax later
	public void printRule()
	{
		System.out.println(left+"->"+Arrays.toString(right)); 
	}
	*/
}
