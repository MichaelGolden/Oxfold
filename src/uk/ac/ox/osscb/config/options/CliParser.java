package uk.ac.ox.osscb.config.options;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class CliParser {
	
	public static void main(String[] args) {
		CliParser cliParser = new CliParser();
		RunOptions parse = cliParser.parse(args);
		System.out.println("got args: " + parse);
	}
	
	
	/**
	 *  &lt;alingment&gt;
	 *  --grammar=ppfold.grammar
	 *  --grammar-params=doc/ppfold.parameters
	 *  --tree-file=doc/newickFiles/TestRNAData16.newick
	 *  --weight=0
	 *  
	 * @param args
	 * @return
	 */
	public RunOptionsPartial parse(String[] args){
        CommandLine line = null;
        
		 // create the parser
	    CommandLineParser parser = new PosixParser();// GnuParser();
	    try {
	        // parse the command line arguments
	        line = parser.parse( buildOptions(), args );
	    }
	    catch( ParseException exp ) {
	        // oops, something went wrong
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );	        
	    }		
	    

	  //  System.err.println("public RunOptionsPartial parse(String[] args) "+OptionsHelper.treeOptName);

	    	    
		String grammar = line.getOptionValue(OptionsHelper.grammarOptName);
		String grammarPath = line.getOptionValue(OptionsHelper.grammarParamsOptName);
		String treePath = line.getOptionValue(OptionsHelper.treeOptName);
		String weightStr = line.getOptionValue(OptionsHelper.weightOptName);
		
		Double weight = OptionsHelper.parseWeight(weightStr, args);
		
	    String[] argList = line.getArgs();
	    String alignmentPath = argList.length > 0 ? argList[0] : null;
	    

		return new RunOptionsPartial(alignmentPath, grammar, grammarPath, treePath, weight);
	}

	private Options buildOptions(){
		Option grammarOpt   = OptionBuilder.withArgName( "grammar definition file" )
                .hasArg()
                .withDescription(  "path to grammar file" )
                .create( OptionsHelper.grammarOptName );
		Option grammarParamsOpt   = OptionBuilder
			.withLongOpt(OptionsHelper.grammarParamsOptName)
			.withDescription(  "path to grammar parameters file" )
			.hasArg()
			.withArgName("grammar_params")
			.create();

		Option evolTreeOpt   = OptionBuilder.withArgName( "evolutionary tree definition file" )
                .hasArg()
                .withDescription( "path to evolutionary tree definition" )
                .create( OptionsHelper.treeOptName );
		
		Option weightOpt   = OptionBuilder.withArgName( OptionsHelper.weightOptName )
                .hasArg()
                .withDescription( "value of weight" )
                .create( OptionsHelper.weightOptName );
		
		Options options = new Options();
		options.addOption(grammarOpt);
		options.addOption(grammarParamsOpt);
		options.addOption(evolTreeOpt);
		options.addOption(weightOpt);
		
		return options;
	}
}
