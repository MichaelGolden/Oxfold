package uk.ac.ox.osscb.grammar;

/**
 * Type of rule.
 *  
 * @author Vladimir
 */
public enum RuleType {
	
	/**
	 *  A->. (aka U->.)
	 */
	RULE1,
	
	/**
	 * A -> (B) (aka U->(V))
	 */
	RULE2,
	
	/**
	 * A->BC (aka U->VW)
	 */
	RULE3
}
