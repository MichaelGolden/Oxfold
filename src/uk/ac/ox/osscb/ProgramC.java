package uk.ac.ox.osscb;

import org.apache.commons.lang3.ClassUtils;

import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

/**
 * 
 * Co-transcriptional folding. Does not seem to work well.
 * 
 * @author Pierre
 */
public class ProgramC {

	public static void main(String[] args) {
		new ProgramC().run(args);
	}

	private void run(String[] args) {
		if(4 != args.length && 5 != args.length){
			printUsage();
			return;
		}

		long start = System.nanoTime();// use nano to measure elapsed time (rather than wall-clock).
		
		double weight = 0;
		String weightStr = null;
				
		if (5 == args.length) {
			weightStr = args[4];
			try {
				weight = Double.valueOf(weightStr.trim());
			} catch (NumberFormatException ex) {
				System.out.println(String.format("Could not parse weight from: '%s'", weightStr));
				printUsage();
				return;
			}
		} else if(4 == args.length) {
			weight = Constants.DefaultWeightParam;
		}

		String alignmentsFile = args[0];
		String grammarFile = args[1];
		String paramsFile = args[2];
		String treeFile = args[3];
		System.out.println(String.format("%s- Running OXFOLD (%s) with Arguments:%n" +
				"\talingment file:\t%s%n\tgrammar definition:\t%s%n\tparameters file:\t%s%n\ttree file:\t%s %n\tweight: %f",
				Util.dateTime(), ClassUtils.getShortClassName(getClass()),
				alignmentsFile, grammarFile, paramsFile, treeFile, weight));
		int alignLen = DefaultAlignmentParser.calculateAlignmentLength(Util.assertCanReadFile(alignmentsFile));
		System.out.println(String.format("\tseq len: %d", alignLen));

		/*
		if(paramsLog.isInfoEnabled()){
			paramsLog.info();
			paramsLog.info(String.format("\talign len: %d", alignLen));
		}*/
		
		new CotranscriptionalFold().newFoldEvolutionary(alignmentsFile,grammarFile,paramsFile,treeFile,weight);
				
		long end = System.nanoTime();

		System.out.println(String.format("%n%1$s- Finished folding. %2$s passed (%3$.2f seconds)",
				Util.dateTime(), Util.spentTimeFromNanoSec(end - start), (end-start)/(double)(1e9)));
	}

	private static void printUsage() {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("Usage:%n\t").
			append(String.format("java %s alignmentsFile grammarFile paramsFile treeFile [weight]%n", Program.class.getName())).
			append("(weight is optional)");
		
		System.out.println(sb);
	}

}
