package uk.ac.ox.osscb.rules;

import java.util.List;

public interface RulesType {

	List<Rule> getRules();
	
	double getAllRulesProbabilities();
}
