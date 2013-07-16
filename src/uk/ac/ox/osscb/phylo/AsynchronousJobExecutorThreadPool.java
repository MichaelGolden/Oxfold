package uk.ac.ox.osscb.phylo;


import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Asynchronous job execution, thread pool class. Jobs are executed in a fixed
 * thread pool.
 * 
 * @author M.Vaerum
 */

public class AsynchronousJobExecutorThreadPool extends AsynchronousJobExecutor {
	private ExecutorService threadPool;
	private int threadCount;

	public AsynchronousJobExecutorThreadPool(int threadCount) {
		this.threadCount = threadCount;
		threadPool = Executors.newFixedThreadPool(threadCount);
	}

	@Override
	public void startExecution(final PhyloJob job, final JobListener listener) {
		threadPool.execute(new Runnable() {
			public void run() {
				if (job.isType() == false) {
					double[][] res = PhyloCalc.calcSingleColumn(job);
					listener.jobFinished(res);
				} else {
					double[][] res = PhyloCalc.calcDoubleColumn(job);
					listener.jobFinished(res);
				}
			}
		});
	}

	@Override
	public String getDescription() {
		return "Local thread execution. Bounded by " + threadCount + " threads";
	}

	@Override
	public String getId() {
		return this.getClass().getName() + threadCount;
	}

	public void shutDown(){
		this.threadPool.shutdownNow(); 
	}

	public boolean isTerminated() {
		return this.threadPool.isTerminated();
	}
	
	// @Override
	// public AsynchronousJobExecutor createJobExecutor(AlgoParameters param) {
	// //we don't have state
	// return this;
	// }

}
