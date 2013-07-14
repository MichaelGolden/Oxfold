package uk.ac.ox.osscb;


import java.util.HashMap;

import uk.ac.ox.osscb.domain.NucleotideProbsConverter;
import uk.ac.ox.osscb.domain.NucleotideProbsDouble;
import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;

/**
 * completely rewritten implementation of NucleotideBasePairingCalculator
 * suggested syntax for baseUnPairedProbs, basePairedProbs: sequence of "bases" is -,A,C,G,U 
 * 
 * @author Vladimir, lepuslapis
 *
 */
public class NucleotideBasePairingProbsCalculator {

	public NucleotideProbsPrecise calculate(int[][] alignment, NucleotideProbsDouble baseProbabilities){
	
		// we are converting to PointRes once in order to avoid convertion done iteratively in a cycle.
		NucleotideProbsPrecise baseProbsPrecise = new NucleotideProbsConverter().toPrecise(baseProbabilities);
		baseProbabilities = null;// to ensure we're not using it in this method.

		int dim = alignment[0].length;
		int seqnum = alignment.length;
		NucleotideProbsPrecise alignmentProbs = new NucleotideProbsPrecise(dim);
		for (int j = 0; j < dim; j++){
			//n[j,0] = product(sbp(a[:,j])); 
			PointRes tmp = PointRes.ONE;
			for (int k=0; k<seqnum; k++) {
				// tmp = tmp.multiply(PointRes.valueOf(baseProbabilities.getUnpairingProbability(alignment[k][j])));
				tmp = tmp.multiply(baseProbsPrecise.getUnpairingProbability(alignment[k][j]));
			}
			alignmentProbs.setUnpairingProbability(j,tmp);
		}
		for(int j = 0; j<dim; j++){
			for(int k = j + 1; k < dim; k++){
				//n[j,k] = product(pbp(a[:,j],a[:,k]));
				PointRes tmp = PointRes.ONE;
				for (int l=0; l<seqnum; l++) {
					//tmp = tmp.multiply(PointRes.valueOf(baseProbabilities.getPairingProbability(alignment[l][j], alignment[l][k])));
					tmp = tmp.multiply(baseProbsPrecise.getPairingProbability(alignment[l][j], alignment[l][k]));
				}
				alignmentProbs.setPairingProbability(j, k, tmp);
			}
		}
		return alignmentProbs;		
	}
	
	public NucleotideProbsPrecise calculate(String[] align, EvolutionaryParameters parameters){
		if(null == align){
			throw new IllegalArgumentException("align array cannot be null");
		}
		if(align.length < 1){
			throw new IllegalArgumentException("align array cannot be empty");
		}


		NucleotideProbsPrecise res = new NucleotideProbsPrecise(align[0].length());
		
		HashMap<String, String[]> sAlphaSynonyms = parameters.getSAlphabet().getSynonyms();
		String[] standardNames = parameters.getSAlphabet().getStandardNames();
		Qmatrix sQmtx = parameters.getSQmatrix();
		
		HashMap<String,Double> sIndices = new HashMap<String,Double>();
		for (int j = 0; j<standardNames.length;j++) {
			sIndices.put(standardNames[j], sQmtx.getPrior()[j]);
		}
			
		
		///////////////////////////////////////////////////
		// unpairing:		
		for(int colIdx = 0; colIdx < align[0].length(); colIdx++){
			
			PointRes unpairingProb = PointRes.valueOf(1);
			
			for(int rowIdx = 0; rowIdx < align.length; rowIdx++){
				String row = align[rowIdx];
				String charStr = row.substring(colIdx, colIdx+1);
				
				String[] strings = sAlphaSynonyms.get(charStr);
				
				double elementPrior = 0;
				for(String syn : strings){
					elementPrior += sIndices.get(syn);
				}
								
				unpairingProb = unpairingProb.multiply(PointRes.valueOf(elementPrior));
			}
			
			res.setUnpairingProbability(colIdx, unpairingProb);
		}
		

		///////////////////////////////////////////////////
		// pairing:
		Alphabet pAb = parameters.getPAlphabet();
		String[] pStdNames = pAb.getStandardNames();
		HashMap<String, String[]> pAbSynonyms = pAb.getSynonyms();
		Qmatrix pQmtx = parameters.getPQmatrix();
		
		HashMap<String, Double> pIndices = new HashMap<String, Double>();
		for (int j = 0; j<pStdNames.length; j++) {
			pIndices.put(pStdNames[j], pQmtx.getPrior()[j]);
		}
		
		int seqnum = align.length;
		for(int lColIdx = 0; lColIdx < align[0].length()-1; lColIdx++){
			
			for(int rColIdx = lColIdx+1; rColIdx < align[0].length(); rColIdx++){
				PointRes pairingProb = PointRes.valueOf(1);
				int unsure = 0;
				for(int rowIdx = 0; rowIdx < align.length; rowIdx++){
					String row = align[rowIdx];
					String charDuplex = row.substring(lColIdx, lColIdx+1) + 
							row.substring(rColIdx, rColIdx+1);
					
					String[] synStrings = pAbSynonyms.get(charDuplex);
					if (synStrings.length > 1) {
						unsure++;
					}
					
					double elementPrior = 0;
					for(String syn : synStrings){
						elementPrior += pIndices.get(syn);
					}
					
					pairingProb = pairingProb.multiply(PointRes.valueOf(elementPrior));
				}
				if (unsure>0.5*seqnum) {
					res.setPairingProbability(lColIdx, rColIdx, PointRes.ZERO);
				} else {
					res.setPairingProbability(lColIdx, rColIdx, pairingProb);
				}
			}
		}
			
		return res;
	}
}
