package uk.ac.ox.osscb.vienna;

import uk.ac.ox.osscb.parser.DefaultAlignmentParser;

public class Fasta2StockholmConvertingViennaRnaRunner implements ViennaRunner{
	
	private static final String stockFileName = "stock.data.file";
	
	private OsCommandBasedViennaRnaRunner osCommandBasedViennaRnaRunner;

	public Fasta2StockholmConvertingViennaRnaRunner(String path2Vienna, String fastaFile) {
		super();
		
		SimpleFasta2StockholmConverter converter = new SimpleFasta2StockholmConverter(fastaFile);
		
		converter.convert(stockFileName);
		
		// hardcore way of dependency injection
		int seqLen = DefaultAlignmentParser.calculateAlignmentLength(fastaFile);

		// hardcore way of dependency injection
		this.osCommandBasedViennaRnaRunner = new 
				OsCommandBasedViennaRnaRunner(path2Vienna, stockFileName, seqLen);
	}

	@Override
	public double[][] getViennaProbabilities(String constrainingStructure) {
		return this.osCommandBasedViennaRnaRunner.getViennaProbabilities(constrainingStructure);
	}
}
