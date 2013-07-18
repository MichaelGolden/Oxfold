package uk.ac.ox.osscb;



/**
 * Console entry point
 * 
 * @author Vladimir
 */
public class Program extends BaseProgram {
	
	public static void main(String[] args) {
		new Program().run(args);
	}

	protected void runIternal(double weight, String treeFile, boolean haveTree,
			String alignmentsFile, String grammarFile, String paramsFile) {
		
		if (haveTree) {
			new KineticFold2().foldEvolutionary(alignmentsFile,grammarFile,paramsFile,treeFile,weight, 1.0);
		} else {
			new KineticFold2().fold(alignmentsFile, grammarFile, paramsFile, weight);
		}
	}

	public boolean isArgsNumberCorrect(String[] args) {
		return args.length >=1 && args.length <= 5; 
	}

	public String getUsageMsg() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("Usage:%n\t").
			append(String.format("java %s alignmentsFile grammarFile paramsFile [treeFile] [weight]%n", Program.class.getName())).
			append("(treeFile and weight are optional; if not treeFile is provided, a non-evolutionary algorithm is run)");
		
		return sb.toString();
	}
}
