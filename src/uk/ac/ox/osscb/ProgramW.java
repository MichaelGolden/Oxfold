package uk.ac.ox.osscb;

import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class ProgramW {

	public static void main(String[] args) {
		String grammarFile = args[1];
		String paramsFile = args[2];
		String treeFile = args[3];
		String alignmentsFile = args[0];
		
		long start = System.nanoTime();// use nano to measure elapsed time (rather than wall-clock).
		
		ProgramOutput.outMsg2(String.format("%s- Running OXFOLD with Arguments:%n" +
				"\talingment file:\t%s%n\tgrammar definition:\t%s%n\tparameters file:\t%s%n\ttree file:\t%s",
				Util.dateTime(), alignmentsFile, grammarFile, paramsFile, treeFile));
		int alignLen = DefaultAlignmentParser.calculateAlignmentLength(Util.assertCanReadFile(alignmentsFile));
		ProgramOutput.outMsg2(String.format("\tseq len: %d", alignLen));

		long end = System.nanoTime();
		
		new WeightedMEA().foldEvolutionary(alignmentsFile, grammarFile, paramsFile, treeFile);
		
		ProgramOutput.outMsg2(String.format("%n%1$s- Finished folding. %2$s passed (%3$.2f seconds)",
				Util.dateTime(), Util.spentTimeFromNanoSec(end - start), (end-start)/(double)(1e9)));
	}
}
