package uk.ac.ox.osscb.config.options;

public interface RunOptions {

	public abstract String getAlignmentPath();

	public abstract String getGrammarPath();

	public abstract String getGrammarParamsPath();

	public abstract String getTreeDefinitionPath();

	public abstract double getWeight();
	
	public boolean hasTree();

}
