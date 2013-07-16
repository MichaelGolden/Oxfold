package uk.ac.ox.osscb.phylo;


import java.util.List;

/**
 * Asynchronous job execution, blocking class. Jobs are executed sequentially.
 * 
 * @author M.Vaerum
 */

public class AsynchronousJobExecutorBlocking extends AsynchronousJobExecutor {
	@Override
	public String getDescription() {
		return "Serial (blocking) execution";
	}
	
	@Override
	public void startExecution(PhyloJob job, JobListener listener) {
		if (job.isType() == false) {
			double[][] res = PhyloCalc.calcSingleColumn(job);
			listener.jobFinished(res);
		} else {
			double[][] res = PhyloCalc.calcDoubleColumn(job);
			listener.jobFinished(res);
		}
	}

	@Override
	public String getId() {
		return this.getClass().getName();
	}

	public void shutDown(){
		//Nothing happens here as we only have one thread. 
	}
	
	public boolean isTerminated() {
		return false;
	}
	
	// @Override
	// public AsynchronousJobExecutor createJobExecutor(AlgoParameters param) {
	// //We have no state
	// return this;
	// }

}
