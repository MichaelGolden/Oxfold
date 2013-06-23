package uk.ac.ox.osscb;

import uk.ac.ox.osscb.inoutside.PPOutputDouble;
import uk.ac.ox.osscb.parser.DefaultAlignmentParser;
import uk.ac.ox.osscb.vienna.ViennaOutputCalculator;

public class ViennaKineticFold {

	//private final Logger log = LoggerFactory.getLogger(KineticFold.class);
	
	public void fold(String alignmentFile){
		Util.assertCanReadFile(alignmentFile);
		
		new DefaultAlignmentParser();
		int length = DefaultAlignmentParser.calculateAlignmentLength(alignmentFile);
		
		// by default is initialised with zeros automatically
		int[] structure = new int[length];
		for(int posIdx = 0; posIdx < structure.length; posIdx++){
			structure[posIdx] = Constants.UnpairedBaseIdx;
		}
				
		boolean exitBecauseOfDiff = false;
		for(int iterSoFar = 0; iterSoFar < Constants.MaxIterations; iterSoFar++){
			PPOutputDouble ppProbs = new ViennaOutputCalculator().calculate(alignmentFile, structure);
			
			if (ppProbs.getDiff()>Constants.IterationCutOffDouble) {
				structure = new StructureUtils().makeNewStructure(structure, ppProbs.getHelix());
				dumpCurrentOutput(ppProbs);
			} else {
				exitBecauseOfDiff = true;
				dumpCurrentOutput(ppProbs);
				break;
			}
			dumpStructure(structure);
		}


		dumpExitReason(exitBecauseOfDiff);
		
		OutputGenerator outputGenerator = new LoggingOutputGenerator();
		outputGenerator.generate(structure);

	}
	private void dumpExitReason(boolean exitBecauseOfDiff) {
		ProgramOutput.outMsg(String.format("\texiting because of %s", 
				exitBecauseOfDiff ? "difference under threshold." : "number of iterations."));
	}

	private void dumpCurrentOutput(PPOutputDouble ppProbs) {
		String msg = String.format("lIdx: %d, rIdx: %d\tHelix Length: %d\tdiff: %g"
				, ppProbs.getLeftIdx()
				, ppProbs.getRightIdx()
				, ppProbs.gethelixLength()
				, ppProbs.getDiff()
				);
		ProgramOutput.outMsg(msg);
	}

	private void dumpStructure(int[] structure) {
		ProgramOutput.outMsg(new LoggingOutputGenerator().dumpStructure(structure));
	}

}
