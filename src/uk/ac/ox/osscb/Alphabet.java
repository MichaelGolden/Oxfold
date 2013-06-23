package uk.ac.ox.osscb;

import java.util.HashMap;

public class Alphabet {
		
		private String[] names;
		private String[] standardNames;
		private HashMap<String,String[]> synonyms;
				
		public Alphabet(Letter[] letters) {
			int length = letters.length;
			int countStandard = 0;
			this.names = new String[length];
			synonyms = new HashMap<String, String[]>();
			for (int j = 0; j<length; j++) {
				this.names[j] = letters[j].getName();
				this.synonyms.put(letters[j].getName(), letters[j].getSynonyms());
				if (letters[j].getIsStandard()) {
					countStandard++;
				}
			}
			String[] tmp = new String[countStandard];
			countStandard = 0;
			for (int j = 0; j<length; j++) {
				if (letters[j].getIsStandard()) {
					tmp[countStandard] = letters[j].getName();
					countStandard++;
				}	
			}
			this.standardNames = tmp;
		}
		
		public String[] getNames() {
			return this.names;
		}
		
		public String[] getStandardNames() {
			return this.standardNames;
		}
		
		public HashMap<String,String[]> getSynonyms() {
			return this.synonyms;
		}
}
