package uk.ac.ox.osscb;

public class EvolutionaryParameters {
		
		/**
		 * sAlphabet - alphabet of single nucleotides
		 * sQmatrix - associated Qmatrix object
		 * pAlphabet - alphabet of paired nucleotides
		 * pQmatrix - associated Qmatrix object
		 */
		private Alphabet sAlphabet;
		private Alphabet pAlphabet;
		private Qmatrix sQmatrix;
		private Qmatrix pQmatrix;
		
		public EvolutionaryParameters(Alphabet sAlphabet, Qmatrix sQmatrix, Alphabet pAlphabet, Qmatrix pQmatrix) {
			this.sAlphabet=sAlphabet;
			this.pAlphabet=pAlphabet;
			this.sQmatrix=sQmatrix;
			this.pQmatrix=pQmatrix;
		}
		
		public Alphabet getSAlphabet() {
			return this.sAlphabet;
		}
		
		public Alphabet getPAlphabet() {
			return this.pAlphabet;
		}
		
		public Qmatrix getSQmatrix() {
			return this.sQmatrix;
		}
		
		public Qmatrix getPQmatrix() {
			return this.pQmatrix;
		}
}
