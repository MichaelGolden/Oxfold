package uk.ac.ox.osscb;

public interface OutputGenerator {
	
	void generate(int[] structure);

	void generateFinal(int[] structure);
	void generate(Structure structure);

	void generateFinal(Structure structure);

}
