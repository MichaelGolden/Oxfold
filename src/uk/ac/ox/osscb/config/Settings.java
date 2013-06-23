package uk.ac.ox.osscb.config;

import org.apache.commons.lang3.ClassUtils;


/**
 * TODO: should this class validate property values or factory?
 * @author Vladimir
 *
 */
public class Settings {
	
	private String viennaPath;
	
	private static Settings instance;

	Settings() {
		super();
	}
	
	static void set(Settings instance){
		Settings.instance = instance;
	}
	
	public static Settings get(){
		if(null == Settings.instance)
			throw new IllegalStateException(String.format(
					"Settings have not been initialised. Forgot to call 'init' method of the %s?",
					ClassUtils.getShortClassName(SettingsFactory.class)));
		return instance;
	}

	void setViennaPath(String viennaPath) {
		// Util.assertCanReadFile(viennaPath);
		this.viennaPath = viennaPath;
	}
	
	public String getViennaPath() {
		return viennaPath;
	}
}
