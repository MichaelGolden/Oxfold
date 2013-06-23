package uk.ac.ox.osscb.rules;

import java.math.BigDecimal;

public interface Rule {

	BigDecimal getRuleProbability();
	
	Character left();
	
	Character[] right();

}
