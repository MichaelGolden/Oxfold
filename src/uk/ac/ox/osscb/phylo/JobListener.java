package uk.ac.ox.osscb.phylo;


import java.util.List;

/**
 * Job listener interface.
 * 
 * @author M.Vaerum
 */
public interface JobListener {
	/**
	 * Used in the case of PhyloJob
	 */
	void jobFinished(double[][] result);

}
