package uk.ac.ox.osscb.config;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class SettingsFactory {
	
	public static final String propViennaPath = "path.vienna";
	
	// private static Settings settings;
	private static String settingsPath;
	
	/**
	 * @deprecated Use {@link Settings#get()} instead
	 * @return
	 */
	@Deprecated
	public static Settings getSettings() {
		return Settings.get();
//		if(null == settings){
//			throw new IllegalStateException("Settings have not been initialised. Forgot to call 'init' method?");
//		}
//		return settings;
	}

	public static Settings init(String propertiesPath){
		
		SettingsFactory.settingsPath = propertiesPath;
		
		Properties properties2 = readProperties(propertiesPath);
		
		Settings settings = new Settings();
		settings.setViennaPath(properties2.getProperty(propViennaPath));
		
//		SettingsFactory.settings = settings;
		Settings.set(settings);
		
		return getSettings();
	}

	/**
	 * @param propPath relative to classpath root. 
	 * Loaded by {@link ClassLoader#getResourceAsStream(String)}
	 * @return
	 */
	private static Properties readProperties(String propPath) {
		
		InputStream res = SettingsFactory.class.getClassLoader().getResourceAsStream(propPath);
		if(null == res){
			throw new IllegalArgumentException(String.format("Could not find resource by path: '%s'",
					propPath));
		}
		
		try{
			
//			FileInputStream fileInputStream = null;
//			fileInputStream = new FileInputStream(propertiesPath);
//			properties.load(fileInputStream);
			
			Properties properties = new Properties();
			properties.load(res);			
			return properties;
			
		} catch (IOException ex) {
			throw new RuntimeException(String.format("Could not read properties from file: '%s'. Error: %s",
					propPath, ex.getMessage()), ex);
		}finally{
			try {
				res.close();
				//fileInputStream.close();
			} catch (IOException ex) {
				throw new RuntimeException(String.format("Could not read properties from file: '%s'. Error: %s", 
						propPath, ex.getMessage()), ex);
			}
		}
	}
}
