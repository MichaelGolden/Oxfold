package uk.ac.ox.osscb.domain;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.MathContext;

import uk.ac.ox.osscb.PointRes;
import uk.ac.ox.osscb.util.ProbabilityValueValidator;

public class NucleotideProbsPrecise extends SquareMatrixPlusVector<PointRes> {
	
	private MathContext mathCtx = null;//new MathContext(5);


	public NucleotideProbsPrecise(int dim) {
		super(dim, PointRes.ZERO);
	}

	@Override
	protected PointRes getInitialNumber() {
		return PointRes.ZERO;
	}

	@Override
	public PointRes setUnpairingProbability(int i, PointRes prob) {
		ProbabilityValueValidator.validateP(prob);
		if(null != this.mathCtx){
			prob = prob.round(this.mathCtx);
		}
		return super.setUnpairingProbability(i, prob);
	}

	@Override
	public PointRes setPairingProbability(int i, int j, PointRes prob) {
		ProbabilityValueValidator.validateP(prob);
		if(null != this.mathCtx){
			prob = prob.round(this.mathCtx);
		}
		return super.setPairingProbability(i, j, prob);
	}
	
	public void writeEvolutionaryProbs(File outFile)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
			for(int i = 0 ; i < this.getDim() ; i++)
			{
				for(int j = 0 ; j <  this.getDim() ; j++)
				{
					writer.write(this.getPairingProbability(i, j)+"\t");
				}
				writer.newLine();
			}

			writer.close();
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		}
	}
	
	
	
	/* 
	public PointRes getUnpairedColumnProduct(int j){
	checkDimensions(j, -1);

	PointRes product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
}

public PointRes getPairedColumnProduct(int j, int k){
	checkDimensions(j, k);
	PointRes product = this.probs[0][j];
	for(int rowIdx = 1; rowIdx < this.dim; rowIdx++){
		product = product.multiply(this.probs[rowIdx][j]);
	}
	return product;
	}
*/
	
}
