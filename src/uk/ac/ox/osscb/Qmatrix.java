package uk.ac.ox.osscb;

public class Qmatrix {

		/**
		 * R - diagonalisation matrix
		 * D - diagonal matrix of eigenvalues
		 * RI - inverse of diagonalisation matrix
		 * prior - prior distribution; stationary distribution of rate matrix, I think.
		 */
		private double[][] rmatrix;
		private double[][] dmatrix;
		private double[][] rImatrix;
		private double[] prior;
		private int size;
	
		public Qmatrix(double[][] rmatrix,double[][] dmatrix,double[][] rImatrix, double[] prior) {
			this.rmatrix = rmatrix;
			this.dmatrix = dmatrix;
			this.rImatrix = rImatrix;
			this.prior = prior;
			this.size=dmatrix.length;
		}
		
		public int getSize() {
			return this.size;
		}
		
		public double[][] getRmatrix() {
			return this.rmatrix;
		}
		
		public double[][] getDmatrix() {
			return this.dmatrix;
		}
		
		public double[][] getRImatrix() {
			return this.rImatrix;
		}
		
		public double[] getPrior() {
			return this.prior;
		}
		
}
