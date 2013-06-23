package uk.ac.ox.osscb.config.options;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import uk.ac.ox.osscb.Constants;

public class RunOptionsFactory {
	
	private List<RunOptionsPartial> options = new LinkedList<RunOptionsPartial>();

	public RunOptionsFactory(String[] args) {
		super();
		RunOptionsPartial cliOpts = new CliParser().parse(args);
		RunOptionsPartial fileOpts = new OptionsFileParser(Constants.MAIN_PROPERTIES_FILE).parse();
		this.options.add(fileOpts);
		this.options.add(cliOpts);
	}

	public void addOptionsPartial(RunOptionsPartial optionsPartial){
		this.options.add(optionsPartial);
	}

	public RunOptions getOpts(){
		String alignmentPath = null;
		String grammarPath = null;
		String grammarParamsPath = null;
		String treePath = null;
		Double weight = null;
		
		for(RunOptionsPartial optPartial : this.options){
			alignmentPath = acceptNewOptValue(alignmentPath, optPartial.getAlignmentPath());
			grammarPath = acceptNewOptValue(grammarPath, optPartial.getGrammarPath());
			grammarParamsPath = acceptNewOptValue(grammarParamsPath, optPartial.getGrammarParamsPath());
			treePath = acceptNewOptValue(treePath, optPartial.getTreeDefinitionPath());
			weight = acceptNewOptValue(weight, optPartial.getWeightNull());
		}
		if(null == weight){
			weight = Constants.DefaultWeightParam;
		}
		
		return new RunOptionsFixed(alignmentPath, grammarPath, grammarParamsPath, treePath, weight);
	}
	
	private static String acceptNewOptValue(String oldOptValue, String newOptValue){
		if(!StringUtils.isBlank(newOptValue)){
			oldOptValue = newOptValue;
		}
		return oldOptValue;
	}
	
	private static Double acceptNewOptValue(Double oldOptValue, Double newOptValue){
		if(null != newOptValue){
			oldOptValue = newOptValue;
		}
		return oldOptValue;
	}
	
	public class RunOptionsFixed implements RunOptions{

		private String alignmentPath;
		private String grammarPath;
		private String grammarParamsPath;
		private String treePath;
		private double weight;
		
		public RunOptionsFixed(String alignmentPath, String grammarPath,
				String grammarParamsPath, String treePath, double weight) {
			super();
			this.alignmentPath = alignmentPath;
			this.grammarPath = grammarPath;
			this.grammarParamsPath = grammarParamsPath;
			this.treePath = treePath;
			this.weight = weight;
		}

		@Override
		public String getAlignmentPath() {
			return this.alignmentPath;
		}

		@Override
		public String getGrammarPath() {
			return this.grammarPath;
		}

		@Override
		public String getGrammarParamsPath() {
			return this.grammarParamsPath;
		}

		@Override
		public String getTreeDefinitionPath() {
			return this.treePath;
		}

		@Override
		public double getWeight() {
			return this.weight;
		}

		@Override
		public boolean hasTree() {
			return !StringUtils.isBlank(getTreeDefinitionPath());
		}
	}
}
