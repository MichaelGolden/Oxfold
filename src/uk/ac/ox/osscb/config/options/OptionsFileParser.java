package uk.ac.ox.osscb.config.options;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class OptionsFileParser {

	private Properties properties = new Properties();

	public OptionsFileParser(String fileName) {
		super();
		
		InputStream is = null;
		try {
			is = new FileInputStream(fileName);
		} catch (FileNotFoundException e) {
		}
		if(null == is){
			is = getClass().getClassLoader().getResourceAsStream(fileName);
		}
		if(null != is){
			try {
				this.properties.load(is);
			} catch (IOException e1) {
			}
			try {
				is.close();
			} catch (IOException e) {
			}
		}
	}

	public RunOptionsPartial parse(){
		// String alignmentPath = properties.getProperty(OptionConstants.) ;
		String grammarPath = this.properties.getProperty(OptionsHelper.grammarOptName);
		String grammarParamsPath = this.properties.getProperty(OptionsHelper.grammarParamsOptName);
		String treeDefinitionPath = this.properties.getProperty(OptionsHelper.treeOptName);
		String weightStr = this.properties.getProperty(OptionsHelper.weightOptName);
		Double Weight = OptionsHelper.parseWeight(weightStr, null);
		
		return new RunOptionsPartial(null, grammarPath, grammarParamsPath, treeDefinitionPath, Weight); 
	}
}
