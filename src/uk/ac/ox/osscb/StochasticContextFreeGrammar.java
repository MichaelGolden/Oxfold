package uk.ac.ox.osscb;


import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import uk.ac.ox.osscb.grammar.Grammar;
import uk.ac.ox.osscb.grammar.RuleType;

/**
 * Date: 16 July 2012
 * @author lmath@cs.ubc.ca
 * Context Free Grammar class
 *
 */

public class StochasticContextFreeGrammar implements Grammar {
	
	private char[] nonterminals;
	private char[] terminals;
	
	//Contains Production Rules of type 1: A-> . 
	private List<ProductionRule> class1 = new LinkedList<ProductionRule>();
	//Contains Production Rules of type 2: A-> (B)
	private List<ProductionRule> class2 = new LinkedList<ProductionRule>();
	//Contains Production Rules of type 2: A-> BC
	private List<ProductionRule> class3 = new LinkedList<ProductionRule>();
	

	/**
	 * Constructor
	 * Parameters must be given when object is created as they do not change.
	 */
	public StochasticContextFreeGrammar(char[] nonterminals, char[] terminals,
			List<ProductionRule> rules)	{
		
		this.nonterminals = nonterminals;
		this.terminals = terminals;
			
		setRuleClasses(rules);
	}
	
	/*
	 * Helper method to sort rules into the rule classes
	 */
	private void setRuleClasses(List<ProductionRule> rules)
	{
		for(ProductionRule prule : rules)
		{
			char[] right = prule.getRight();
			char firstchar = right[0];
			char endchar = right[right.length - 1];

			// Set type
			//TODO: add error handling if a rule is not the correct type
			// class 1: A->.
			if (right.length == 1 && isTerminal(firstchar, terminals))
			{
				prule.ruleclass = RuleType.RULE1;
				class1.add(prule);
			}
			//class 2: A -> (B)
			else if (isTerminal(firstchar, terminals) && isTerminal(endchar, terminals) && right.length > 2)
			{
				prule.ruleclass = RuleType.RULE2;
				class2.add(prule);
			}
			//class 3: A->BC
			else 
			{
				prule.ruleclass = RuleType.RULE3;
				class3.add(prule);
			}
			
		}
	}
	
	/**
	 * Getter for the set of nonterminals
	 * @return nonterminals
	 */
	public char[] getNonterminals()
	{
		return nonterminals;
	}
	
	/**
	 * Method to get a vector of ProductionRules for a rule class 1, 2, or 3
	 * @param ruleclass
	 * @return the class of rules for the int ruleclass specified, or null if ruleclass
	 * not 1, 2, or 3
	 */
	public List<ProductionRule> getRules(RuleType ruleclass)
	{
		switch(ruleclass){
			case RULE1:
				return Collections.unmodifiableList(this.class1);
			case RULE2:
				return Collections.unmodifiableList(this.class2);
			case RULE3:
				return Collections.unmodifiableList(this.class3);
			default:
				throw new IllegalArgumentException(String.format("rule type %s is not supported", ruleclass));
		}
	}
	
	/* (non-Javadoc)
	 * @see uk.ac.osscb.Grammar#getProbabilities()
	 */

	/**
	 * Method to return an array of probabilities for the ruleclass specifed
	 * @param ruleclass -- class of rules 1, 2, 3 as defined in ProductionRule.java
	 * @return an array of probabilities for the specified ruleclass
	 */
	public PointRes[] getProbabilities(RuleType ruleclass)
	{
		List<ProductionRule> rulesforRuleclass = getRules(ruleclass);
		PointRes[] probs = new PointRes[rulesforRuleclass.size()];
		
		for(int i=0; i<rulesforRuleclass.size(); i++)
			probs[i] = rulesforRuleclass.get(i).getProbability();
			
		return probs;
	}


	private boolean isTerminal(char ch, char[] terminals)
	{
		for(int i=0; i<terminals.length; i++)
		{
			if (terminals[i] == ch)
			{
				return true;
			}
		}
		return false;
	}
}
