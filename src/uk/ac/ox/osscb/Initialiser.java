package uk.ac.ox.osscb;

import uk.ac.ox.osscb.config.SettingsFactory;

/**
 * Smth that should be done via IoC, but we're lazy for a while
 * 
 * @author Vladimir
 *
 */
public class Initialiser {
	
	public static void init(){
		SettingsFactory.init("live.settings.properties");
	}

}
