package uk.ac.ox.osscb.config.options;

import org.apache.commons.lang3.StringUtils;

public class OptionsHelper {

	public static final String grammarOptName = "grammar";
	public static final String grammarParamsOptName = "grammar-params";
	public static final String treeOptName = "tree";
	public static final String weightOptName = "weight";
	
	public static Double parseWeight(String weightStr, String[] args) {
		Double weight = null;
		if(null != weightStr && weightStr.trim().length() > 0){
			try{
				weight = Double.parseDouble(weightStr);
				
			}catch(NumberFormatException nfe){
				String argsStr = null == args ? "" :
					String.format(", arguments were:%n\t%s", StringUtils.join(args, ";\t"));

				throw new RuntimeException(String.format("Failed to parse weight: '%s'%s",
						weightStr,argsStr), nfe);
			}
		}
		return weight;
	}
}
