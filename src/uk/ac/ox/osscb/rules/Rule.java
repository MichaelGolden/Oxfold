package uk.ac.ox.osscb.rules;

import uk.ac.ox.osscb.PointRes;


public interface Rule {

	PointRes getRuleProbability();
	
	Character left();
	
	Character[] right();

}
