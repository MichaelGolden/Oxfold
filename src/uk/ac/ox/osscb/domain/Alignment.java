package uk.ac.ox.osscb.domain;

/**
 * Matrix of dot-braces alignments.
 * 
 * @author Vladimir
 *
 */
public class Alignment {
	
	private final char[][] alignment;

	private final int seqLength;
	
	public Alignment(char[][] alignment) {
		super();
		
		if(null == alignment)
			throw new IllegalArgumentException("TODO");
		
		if(alignment.length > 0){
			this.seqLength = alignment[0].length;
			for(int seqIdx = 1; seqIdx < alignment.length; seqIdx++){
				if(alignment[seqIdx].length != this.seqLength){
					throw new IllegalArgumentException("Todo");
				}
			}
		}else{
			this.seqLength = 0;
		}
		
		this.alignment = alignment;
	}
	
	public String getAlignmentStr(int seqIdx){
		return new String(this.alignment[seqIdx]);
	}
	
	public char[] getAlignment(int seqIdx){
		return this.alignment[seqIdx];
	}
	
	public int getAlignmentsCount(){
		return this.alignment.length;
	}
	
	public int getSequenceLength(){
		return this.seqLength;
	}
}
