package uk.ac.ox.osscb.analysis;
public class BasePairMetrics {
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corrsponding to the predicted structure returns the sensitivity.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 * @return
	 */
	public static double calculateSensitivity (int [] realPairedSites, int [] predictedPairedSites)
	{
		double count = 0;
		double total = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(realPairedSites[i] != 0)
			{
				total++;
				if(realPairedSites[i] == predictedPairedSites[i])
				{
					count++;
				}
			}
		}
		
		if(total == 0)
		{
			return 0;
		}
		
		return count / total;
	}
	
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corresponding to the predicted structure returns the PPV.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 * @return
	 */
	public static double calculatePPV (int [] realPairedSites, int [] predictedPairedSites)
	{
		double count = 0;
		double total = 0;
		for(int i = 0 ; i < realPairedSites.length ; i++)
		{
			if(predictedPairedSites[i] != 0)
			{
				total++;
				if(predictedPairedSites[i] == realPairedSites[i])
				{
					count++;
				}
			}
		}
		
		if(total == 0)
		{
			return 0;
		}
		
		return count / total;
	}
	
	/**
	 * Given an array of paired sites corresponding to the real structure 
	 * and an array corresponding to the predicted structure returns the F-score.
	 * @param realPairedSites
	 * @param predictedPairedSites
	 * @return
	 */
	public static double calculateFScore (int [] realPairedSites, int [] predictedPairedSites)
	{
		double sensitivity = calculateSensitivity(realPairedSites, predictedPairedSites);
		double ppv = calculatePPV(realPairedSites, predictedPairedSites);
		
		double fscore =  (2 * sensitivity * ppv)/(sensitivity+ppv);
		return Double.isNaN(fscore) ? 0 : fscore;
	}
}
