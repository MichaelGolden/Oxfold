package uk.ac.ox.osscb.domain;

import uk.ac.ox.osscb.PointRes;



public class NucleotideProbsConverter {
	
	public NucleotideProbsPrecise toPrecise(NucleotideProbsDouble probsDouble){
		
		if(null == probsDouble)
			throw new IllegalArgumentException("probsDouble must not be null!");
		
		final int dim = probsDouble.getDim();
		
		NucleotideProbsPrecise probsPrecise = new NucleotideProbsPrecise(dim);
		for(int i = 0; i < dim; i++){
						
			for(int j = 0; j < dim; j++){
				Double pairingProbDbl = probsDouble.getPairingProbability(i, j);
				probsPrecise.setPairingProbability(i, j, PointRes.valueOf(pairingProbDbl));
			}
			
			Double unPairingProbDbl = probsDouble.getUnpairingProbability(i); 
			probsPrecise.setUnpairingProbability(i, PointRes.valueOf(unPairingProbDbl));
			
		}
		
		return probsPrecise;
	}

}
