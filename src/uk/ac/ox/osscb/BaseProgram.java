package uk.ac.ox.osscb;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.ClassUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.MDC;

import uk.ac.ox.osscb.config.options.RunOptions;
import uk.ac.ox.osscb.config.options.RunOptionsFactory;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public abstract class BaseProgram {
	
	protected static final Pattern alignmentSuffix = Pattern.compile(".+?(\\d+)\\.dat$", Pattern.CASE_INSENSITIVE);

	public void run(String[] args){
		if(!isArgsNumberCorrect(args)){
			printUsage();
			return;
		}
		long start = System.nanoTime();// use nano to measure elapsed time (rather than wall-clock).
		
		RunOptions runOpts = new RunOptionsFactory(args).getOpts();
		
		String alignmentsFile = runOpts.getAlignmentPath();
		String grammarFile = runOpts.getGrammarPath();
		String paramsFile = runOpts.getGrammarParamsPath();
		
		String outFileSuffix = resolveOutFileSuffix(runOpts.hasTree(), runOpts.getWeight(), alignmentsFile);
		MDC.put(Constants.RUN_SUFFIX, outFileSuffix);
		
		ProgramOutput.outMsg2(String.format("%s- Running OXFOLD (%s) with Arguments:%n" +
				"\talignment file:\t%s%n\tgrammar definition:\t%s%n\tparameters file:\t%s%n\ttree file:\t%s %n\tweight: %.2f" +
				"%n\trun suffix: %s",
				Util.dateTime(), ClassUtils.getShortClassName(getClass()), 
				alignmentsFile, grammarFile, paramsFile, runOpts.getTreeDefinitionPath(),
						runOpts.getWeight(), outFileSuffix));
		int alignLen = DefaultAlignmentParser.calculateAlignmentLength(Util.assertCanReadFile(alignmentsFile));
		ProgramOutput.outMsg2(String.format("\tseq len: %d", alignLen));

		Initialiser.init();

		runIternal(runOpts.getWeight(), runOpts.getTreeDefinitionPath(), runOpts.hasTree(), alignmentsFile, grammarFile, paramsFile);

		long end = System.nanoTime();

		ProgramOutput.outMsg2(String.format("%n%1$s- Finished folding. %2$s passed (%3$.2f seconds)",
				Util.dateTime(), Util.spentTimeFromNanoSec(end - start), (end-start)/(double)(1e9)));	
	}

	public void run_Old(String[] args){
		if(!isArgsNumberCorrect(args)){
			printUsage();
			return;
		}

		long start = System.nanoTime();// use nano to measure elapsed time (rather than wall-clock).
		
		double weight = 0;
		String treeFile = "";
		boolean haveTree = false;
		
		if (5 == args.length) {
			haveTree = true;
			treeFile = args[3];
			String weightStr = args[4];
			try {
				weight = Double.valueOf(weightStr.trim());
			} catch (NumberFormatException ex) {
				ProgramOutput.outMsg2(String.format("Could not parse weight from: '%s'", weightStr));
				printUsage();
				return;
			}
		} else if(4 == args.length) {
			String weightCandidateStr = args[3];
			try {
				weight = Double.valueOf(weightCandidateStr.trim());
			} catch (NumberFormatException ex) {
				haveTree = true;
				treeFile = weightCandidateStr;
				weight = Constants.DefaultWeightParam;
			}
		} else if (3 == args.length) {
			weight = Constants.DefaultWeightParam;
		}

		String alignmentsFile = args[0];
		String grammarFile = args[1];
		String paramsFile = args[2];
		
		String outFileSuffix = resolveOutFileSuffix(haveTree, weight, alignmentsFile);
		MDC.put(Constants.RUN_SUFFIX, outFileSuffix);
		
		ProgramOutput.outMsg2(String.format("%s- Running OXFOLD (%s) with Arguments:%n" +
				"\talingment file:\t%s%n\tgrammar definition:\t%s%n\tparameters file:\t%s%n\ttree file:\t%s %n\tweight: %f" +
				"%n\trun suffix: %s",
				Util.dateTime(), ClassUtils.getShortClassName(getClass()), 
				alignmentsFile, grammarFile, paramsFile, treeFile,
						weight, outFileSuffix));
		int alignLen = DefaultAlignmentParser.calculateAlignmentLength(Util.assertCanReadFile(alignmentsFile));
		ProgramOutput.outMsg2(String.format("\tseq len: %d", alignLen));

		Initialiser.init();

		runIternal(weight, treeFile, haveTree, alignmentsFile, grammarFile, paramsFile);

		long end = System.nanoTime();

		ProgramOutput.outMsg2(String.format("%n%1$s- Finished folding. %2$s passed (%3$.2f seconds)",
				Util.dateTime(), Util.spentTimeFromNanoSec(end - start), (end-start)/(double)(1e9)));	
	}
	
	/**
	 * @deprecated Use apache commons cli to print usage
	 */
	@Deprecated()
	public void printUsage() {
		ProgramOutput.outMsg2(getUsageMsg());
	}

	public abstract String getUsageMsg();
	
	public abstract boolean isArgsNumberCorrect(String[] args);

	protected abstract void runIternal(double weight, String treeFile, boolean haveTree,
			String alignmentsFile, String grammarFile, String paramsFile);
	
	public static String resolveOutFileSuffix(boolean haveTree, double weight, String alignmentsFile){
		String runType;
		if(!haveTree) {
			if(0== weight){
				runType = "A";
			}else if(weight > 0){
				runType = "B";
			} else{
				throw new IllegalArgumentException(String.format("Weight must be NON-negative. Got: %g", weight));
			}
		} else {
			if(0== weight){
				runType = "C";				
			}else if(weight > 0){
				runType = "D";				
			}else{
				throw new IllegalArgumentException(String.format("Weight must be NON-negative. Got: %g", weight));
			}
		}
		
		Matcher matcher = alignmentSuffix.matcher(alignmentsFile);
		if(matcher.matches()){
			String seqNo = matcher.group(1);
			// a trick to pad the output with zeros. Couldn't make String.format
			// zero-pad the %s conversion.
			seqNo = String.format("%1s", seqNo);
			runType = seqNo + runType;
		}
		
		return runType;
	}

}
