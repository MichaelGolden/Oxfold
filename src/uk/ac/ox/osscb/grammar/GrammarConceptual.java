package uk.ac.ox.osscb.grammar;



import uk.ac.ox.osscb.InsideOutsideProbabilities;
import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;

/**
 * Defines contract for interaction with an abstract 'Grammar'.
 * 
 * @author Vladimir
 */
public interface GrammarConceptual {

	/**
	 * Return: 
	 * a) single double - highest prob;
	 * b) two integers: indices of bases to pair;
	 * 
	 * @param insideProbs
	 * @param outsideProbs
	 * @param pairingProbs
	 * @param distances
	 * @param weight
	 */
	GrammarOutput calculateMaxBasingProbabilities(InsideOutsideProbabilities insideProbs, InsideOutsideProbabilities outsideProbs,
			NucleotideProbsPrecise pairingProbs, PointRes[][] distances, double weight);
	
	
	GrammarOutput  doGrammar(NucleotideProbsPrecise pairingProbs, PointRes[][] distances, double weight);
}
