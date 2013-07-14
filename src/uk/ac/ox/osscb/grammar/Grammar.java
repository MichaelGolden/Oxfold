package uk.ac.ox.osscb.grammar;

import java.util.List;

import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.ProductionRule;
import uk.ac.ox.osscb.grammar.RuleType;


public interface Grammar {

	public PointRes[] getProbabilities(RuleType ruleclass);

	public char[] getNonterminals();
	
	public List<ProductionRule> getRules(RuleType ruleclass);
/*
	public Rule[] getRules1();
	
	public Rule[] getRules2();
	
	public Rule[] getRules3();
	*/
	
}
