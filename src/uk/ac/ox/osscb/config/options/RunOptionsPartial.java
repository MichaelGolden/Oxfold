package uk.ac.ox.osscb.config.options;

import org.apache.commons.lang3.StringUtils;

import uk.ac.ox.osscb.Constants;

/**
 * Option which were used to run the program, provided either
 * via CLI or provided in {@link Constants#MAIN_PROPERTIES_FILE} 
 * @author Vladimir
 *
 */
public class RunOptionsPartial implements RunOptions {
	
	private String alignmentPath;
	private String grammarPath;
	private String grammarParamsPath;
	private String treeDefinitionPath;
	private Double weight;
	
	public RunOptionsPartial(String alignmentPath, String grammarPath,
			String grammarParamsPath, String treeDefinitionPath, Double weight) {
		super();
		this.alignmentPath = alignmentPath;
		this.grammarPath = grammarPath;
		this.grammarParamsPath = grammarParamsPath;
		this.treeDefinitionPath = treeDefinitionPath;
		this.weight = weight;
	}


	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.config.options.RunOptionsI#getAlignmentPath()
	 */
	@Override
	public String getAlignmentPath() {
		return alignmentPath;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.config.options.RunOptionsI#getGrammarPath()
	 */
	@Override
	public String getGrammarPath() {
		return grammarPath;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.config.options.RunOptionsI#getGrammarParamsPath()
	 */
	@Override
	public String getGrammarParamsPath() {
		return grammarParamsPath;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.config.options.RunOptionsI#getTreeDefinitionPath()
	 */
	@Override
	public String getTreeDefinitionPath() {
		return treeDefinitionPath;
	}

	/* (non-Javadoc)
	 * @see uk.ac.ox.osscb.config.options.RunOptionsI#getWeight()
	 */
	@Override
	public double getWeight() {
		if(null != this.weight)
			return weight.doubleValue();
		throw new IllegalStateException("weight has not been initialised");
	}
	
	public Double getWeightNull(){
		return this.weight;
	}

	@Override
	public String toString() {
		return "RunOptionsPartial [alignmentPath=" + alignmentPath
				+ ", grammarPath=" + grammarPath + ", grammarParamsPath="
				+ grammarParamsPath + ", treeDefinitionPath="
				+ treeDefinitionPath + ", weight=" + weight + "]";
	}

	@Override
	public boolean hasTree() {
		return !StringUtils.isBlank(getTreeDefinitionPath());
	}
}
