package uk.ac.ox.osscb.vienna;

import uk.ac.ox.osscb.Constants;
import uk.ac.ox.osscb.HelicesMaker;
import uk.ac.ox.osscb.LoggingOutputGenerator;
import uk.ac.ox.osscb.PossiblePairFinder;
import uk.ac.ox.osscb.config.Settings;
import uk.ac.ox.osscb.config.SettingsFactory;
import uk.ac.ox.osscb.inoutside.Helix;
import uk.ac.ox.osscb.inoutside.PPOutputDouble;

public class ViennaOutputCalculator {

		public PPOutputDouble calculate(String alignmentPath, int[] structure) {
			
			String structureString = new LoggingOutputGenerator().dumpStructure(structure);
			boolean[][] canPair = new PossiblePairFinder().canPair(structure); 
			double[][] pairingProbs = 
					//new OsCommandBasedViennaRnaRunner(SettingsFactory.getSettings().getViennaPath(), new SimpleFasta2StockholmConverter(alignmentPath)).getViennaProbabilities(alignmentPath, structureString);
					new Fasta2StockholmConvertingViennaRnaRunner(Settings.get().getViennaPath(), alignmentPath).getViennaProbabilities(structureString);
			double[][] diffs = getDiffs(pairingProbs);
			
			double maxp = -0.5;// doesn't matter much
			
			int leftIdx = -1; int rightIdx = -1;
			for (int j = 0; j < pairingProbs.length; j++) {
				if (Constants.UnpairedBaseIdx == structure[j]) {
					for (int k = j+1; k<pairingProbs.length; k++) {
						if (pairingProbs[j][k]>maxp) {
							leftIdx = j; rightIdx = k; maxp = pairingProbs[j][k];
						}
					}
				}
			}
			double diff = diffs[leftIdx][rightIdx];
			Helix helix = new HelicesMaker().makeHelix(leftIdx, rightIdx, diffs, canPair);
			return new PPOutputDouble(helix, diff);
		}
		
		private double[] getUnpairingProbs(double[][] pairingProbs) {
			double[] unpairingProbs = new double[pairingProbs.length];
			for (int j = 0; j<unpairingProbs.length; j++) {
				double tmp = 1.0;
				for (int k = 0; k<j; k++) {
					tmp = tmp - pairingProbs[k][j];
				}
				for (int k = j+1; k<pairingProbs.length; k++) {
					tmp = tmp - pairingProbs[j][k];
				}
			}
			return unpairingProbs;
		}
		
		private double[][] getDiffs(double[][] pairingProbs) {
			double[] unpairingProbs = getUnpairingProbs(pairingProbs);
			double[][] diffs = new double[pairingProbs.length][pairingProbs.length];
			for (int j = 0; j<pairingProbs.length; j++) {
				for (int k = j+1; k<pairingProbs.length; k++) {
					diffs[j][k] = pairingProbs[j][k]-(unpairingProbs[j]+unpairingProbs[k])/2;
				}
			}
			return diffs;
		}
}
