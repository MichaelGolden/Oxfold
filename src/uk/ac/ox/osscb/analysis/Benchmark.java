package uk.ac.ox.osscb.analysis;

public class Benchmark {
	public StructureData experimental;
	public StructureData predicted;
	public double fscore;
	public double sensitivity;
	public double ppv;
	public double mountain;
	
	public double param1;
	
	
	private Benchmark(StructureData experimental, StructureData predicted)
	{
		this.experimental = experimental;
		this.predicted = predicted;
	}
	
	public static Benchmark benchmark(StructureData experimental, StructureData predicted)
	{
		double sensitivity = BasePairMetrics.calculateSensitivity(experimental.pairedSites, predicted.pairedSites);
		double ppv = BasePairMetrics.calculatePPV(experimental.pairedSites, predicted.pairedSites);
		double fscore = BasePairMetrics.calculateFScore(experimental.pairedSites, predicted.pairedSites);
		double mountainSim = MountainMetrics.calculateWeightedMountainSimilarity(experimental.pairedSites, predicted.pairedSites);
		
		Benchmark benchmark = new Benchmark(experimental, predicted);
		benchmark.fscore = fscore;
		benchmark.sensitivity = sensitivity;
		benchmark.ppv = ppv;
		benchmark.mountain = mountainSim;
		
		return benchmark;		
	}
	
	public String toString()
	{
		return experimental.file.getName()+","+experimental.pairedSites.length+","+sensitivity+","+ppv+","+fscore+","+mountain;
	}
}
