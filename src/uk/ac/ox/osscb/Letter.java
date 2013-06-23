package uk.ac.ox.osscb;

public class Letter {

		private String name;
		private String[] synonyms;
		private boolean isStandard;
				
		public Letter(String name, String[] synonyms) {
			this.name = name;
			this.synonyms = synonyms;
			if (1 == this.synonyms.length) {
				this.isStandard = true;
			} else {
				this.isStandard = false;
			}
		}
		
		public String getName() {
			return this.name;
		}
		
		public String[] getSynonyms() {
			return this.synonyms;
		}
		
		public boolean getIsStandard() {
			return this.isStandard;
		}
}
