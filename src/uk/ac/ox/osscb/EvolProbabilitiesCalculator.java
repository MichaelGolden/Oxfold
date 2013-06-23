package uk.ac.ox.osscb;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Stack;

import javax.swing.tree.DefaultMutableTreeNode;

import uk.ac.ox.osscb.domain.NucleotideProbsPrecise;

/**
 * Implementation of postorder traversal + calculation of probs for 
 * evolutionary model
 * @author lmath, lepuslapis
 *
 */
public class EvolProbabilitiesCalculator {
		
	public BigDecimal postOrderStackTraversal(EvolutionaryTree tree, Qmatrix qMatrix,
			Alphabet alphabet, double[] prior, String[] leafData, double[][][] transitionMatrices){
	
		Stack<DefaultMutableTreeNode> stack = new Stack<DefaultMutableTreeNode>();

		String[] standardLetters = alphabet.getStandardNames();
		LinkedList<String> standardLettersList = new LinkedList<String>(Arrays.asList(standardLetters));
		HashMap<String, String[]> synonyms = alphabet.getSynonyms();
		
		int totalNodes = tree.nodesInTreeCount(tree.getRoot());
		//matrix of likelihoods
		BigDecimal[][] L = new BigDecimal[totalNodes][standardLetters.length];
		//just zero L
		for(int i=0; i<totalNodes; i++)
			for(int j=0; j<standardLetters.length; j++)
				L[i][j] = BigDecimal.ZERO;
		
		//count the #leaves we've seen so that we can compare to the array of leafData
		int leafCount = 0;
		//traversal number so that we know where to find entries in L
		int traversalIdx = 0;
		
		DefaultMutableTreeNode currNode = tree.getRoot();
		stack.push(currNode);
		while(!stack.empty()){
			boolean hitLeaf = false;
			
			if(currNode.isLeaf() == false){ //currNode is not a leaf
				//if there are leaf nodes below currNode, push this node onto the stack
				//tree.getRoot().getL
				
				//currNode.getFirstLeaf().getLevel() > currNode.getLevel() && (!stack.contains(currNode)|| currNode.equals(tree.getRoot()))
				if(hitLeaf == false){
					currNode = (DefaultMutableTreeNode) currNode.getFirstChild();
					stack.push(currNode);
				}	
				else{
					//iterate over children to get L
					// speeding up the whole thing by calculating fewer matrices
					// current node is not a leaf, so iterate over its children
					Enumeration children = currNode.children();
					int childrenCount = currNode.getChildCount();
					int[] childIndices = new int[childrenCount];
					int count = 0;
					while (children.hasMoreElements()) {
						DefaultMutableTreeNode nextChild = (DefaultMutableTreeNode) children.nextElement();
						EvolutionaryNodeData childData = (EvolutionaryNodeData) nextChild.getUserObject();
						childIndices[count] = childData.getIdx();
						count++;
					}
					for (int ltrIdx1 = 0; ltrIdx1<standardLetters.length; ltrIdx1++) {
						double tmp = 1.0;
						for (int i=0; i<childrenCount; i++) {
							double tmpp = 0.0;
							for (int ltrIdx2 = 0; ltrIdx2<standardLetters.length; ltrIdx2++) {
								double tmpTransProb = transitionMatrices[i][ltrIdx1][ltrIdx2];
								tmpp = tmpp + L[childIndices[i]][ltrIdx2].doubleValue()*tmpTransProb;
							}
							tmp = tmp * tmpp;
						}
						L[traversalIdx][ltrIdx1] = BigDecimal.valueOf(Math.abs(tmp));
					}
					currNode = stack.pop();
					traversalIdx++;
					hitLeaf = false;
				}
				
			}
			
			else{ //currNode is leaf (never added to stack)
				
				hitLeaf = true;
				String observed = leafData[leafCount];
				for(String obs : synonyms.get(observed))
					L[traversalIdx][standardLettersList.indexOf(obs)] = BigDecimal.ONE;
				
				DefaultMutableTreeNode parent = (DefaultMutableTreeNode) currNode.getParent();
				parent.remove(currNode);
				currNode = parent;
				leafCount++;
				traversalIdx++;
			}
			
		}
		BigDecimal likelihood = BigDecimal.ZERO;
		for(int ltrIdx=0; ltrIdx<standardLetters.length; ltrIdx++)
			likelihood = likelihood.add(L[totalNodes-1][ltrIdx].multiply(BigDecimal.valueOf(prior[ltrIdx])));
		//System.out.println(likelihood.doubleValue());
		return likelihood;
		
	}
	
	
	
	
	/**
	 * Calculates and returns the likelihood of the input tree given the 
	 * residues at the leaves. 
	 * @param tree - a evolutionary tree
	 * @param rateMatrix - rate matrix for transitions between different elements
	 * of the alphabet. 
	 * @param alphabet - alphabet of RNA bases with a gap character
	 * @param prior - prior distribution at the root
	 * @return likelihood
	 */
	public BigDecimal postOrderTraversal(EvolutionaryTree tree, Qmatrix qMatrix,
				Alphabet alphabet, double[] prior, String[] leafData, double[][][] transitionMatrices){
		
		String[] standardLetters = alphabet.getStandardNames();
		LinkedList<String> standardLettersList = new LinkedList<String>(Arrays.asList(standardLetters));
		HashMap<String, String[]> synonyms = alphabet.getSynonyms();
		
		int totalNodes = tree.nodesInTreeCount(tree.getRoot());
		//matrix of likelihoods
		BigDecimal[][] L = new BigDecimal[totalNodes][standardLetters.length];
		//just zero L
		for(int i=0; i<totalNodes; i++)
			for(int j=0; j<standardLetters.length; j++)
				L[i][j] = BigDecimal.ZERO;
				
		//count the #leaves we've seen so that we can compare to the array of leafData
		int leafCount = 0;
		//traversal number so that we know where to find entries in L
		int traversalIdx = 0;
		

		
		//go through a post-order traversal and track the current index so that we can update L
		Enumeration postOrderEnum = tree.getRoot().postorderEnumeration();
		while(postOrderEnum.hasMoreElements()){
			DefaultMutableTreeNode nextNode = (DefaultMutableTreeNode) postOrderEnum.nextElement();
			// add the traversalIdx to the object of nextNode
			((EvolutionaryNodeData) nextNode.getUserObject()).setIdx(traversalIdx);
			//	EvolutionaryNodeData data = (EvolutionaryNodeData) nextNode.getUserObject();
			//	data.setIdx(traversalIdx);
			//	nextNode.setUserObject(data);
								
			//if the node is a leaf, then set the residue and all its synonyms
			//to have probability 1
			if (nextNode.isLeaf()){
				String observed = leafData[leafCount];
				for(String obs : synonyms.get(observed)) {
					L[traversalIdx][standardLettersList.indexOf(obs)] = BigDecimal.ONE;
				}
				leafCount++;
			} else {
				// speeding up the whole thing by calculating fewer matrices
				// current node is not a leaf, so iterate over its children
				Enumeration children = nextNode.children();
				int childrenCount = nextNode.getChildCount();
				int[] childIndices = new int[childrenCount];
				int count = 0;
				while (children.hasMoreElements()) {
					DefaultMutableTreeNode nextChild = (DefaultMutableTreeNode) children.nextElement();
					EvolutionaryNodeData childData = (EvolutionaryNodeData) nextChild.getUserObject();
					childIndices[count] = childData.getIdx();
					count++;
				}
				for (int ltrIdx1 = 0; ltrIdx1<standardLetters.length; ltrIdx1++) {
					double tmp = 1.0;
					for (int i=0; i<childrenCount; i++) {
						double tmpp = 0.0;
						for (int ltrIdx2 = 0; ltrIdx2<standardLetters.length; ltrIdx2++) {
							double tmpTransProb = transitionMatrices[childIndices[i]][ltrIdx1][ltrIdx2];
							tmpp = tmpp + L[childIndices[i]][ltrIdx2].doubleValue()*tmpTransProb;
						}
						tmp = tmp * tmpp;
					}
					L[traversalIdx][ltrIdx1] = BigDecimal.valueOf(Math.abs(tmp));
				}
			}	
			traversalIdx++;
		}
		//String mtx = Util.print2DArray(L);
		//System.out.println(mtx);
		
		BigDecimal likelihood = BigDecimal.ZERO;
		for(int ltrIdx=0; ltrIdx<standardLetters.length; ltrIdx++)
			likelihood = likelihood.add(L[totalNodes-1][ltrIdx].multiply(BigDecimal.valueOf(prior[ltrIdx])));
		//System.out.println(likelihood.doubleValue());
		return likelihood;
	}
	
	public double[][][] getMatrices(EvolutionaryTree tree, Qmatrix qMatrix) {
		Enumeration postOrderEnum = tree.getRoot().postorderEnumeration();
		int matrixSize = qMatrix.getSize();
		int traversalIdx = 0;
		int totalNodes = tree.nodesInTreeCount(tree.getRoot());
		double[][][] result = new double[totalNodes][matrixSize][matrixSize];
		while (postOrderEnum.hasMoreElements()) {
			DefaultMutableTreeNode nextNode = (DefaultMutableTreeNode) postOrderEnum.nextElement();
			((EvolutionaryNodeData) nextNode.getUserObject()).setIdx(traversalIdx);
			if (!nextNode.isLeaf()) {
				Enumeration children = nextNode.children();
				while (children.hasMoreElements()) {
					DefaultMutableTreeNode nextChild = (DefaultMutableTreeNode) children.nextElement();
					EvolutionaryNodeData childData = (EvolutionaryNodeData) nextChild.getUserObject();
					double distance = tree.getDistance(nextNode, nextChild);
					result[childData.getIdx()] = getTransitionProbsMatrix(qMatrix, distance);
				}
			}
			traversalIdx++;
		}
		return result;
	}
	
	public NucleotideProbsPrecise getEvolutionaryProbs(EvolutionaryTree tree, EvolutionaryParameters param, String[] alignment){

		NucleotideProbsPrecise probs = new NucleotideProbsPrecise(alignment[0].length());
		
		//take alignment and turn it into a String[] of columns to input as leafData
		String[] cols = convertToCols(alignment);
		//add unpairing probabilities
		Alphabet sAlphabet = param.getSAlphabet();
		//boolean[] highlyConservedCols = findConservedCols(cols,sAlphabet);
		Qmatrix sQmatrix = param.getSQmatrix();
		double[] sPrior = sQmatrix.getPrior();
		//get the transition matrices first, before calculating them over and over again
		double[][][] sMatrices = getMatrices(tree, sQmatrix);
		for(int colIdx=0; colIdx<cols.length; colIdx++){
			String[] leafData = toLeafData(cols[colIdx]);
			BigDecimal likelihood = postOrderTraversal(tree, sQmatrix, sAlphabet, sPrior, leafData,sMatrices);
			probs.setUnpairingProbability(colIdx, likelihood);
		}
		
		//add pairing probabilities
		Alphabet pAlphabet = param.getPAlphabet();
		Qmatrix pQmatrix = param.getPQmatrix();
		double[] pPrior = pQmatrix.getPrior();
		String[] pStdNames = pAlphabet.getStandardNames();
		HashMap<String, Double> pIndices = new HashMap<String, Double>();
		for (int j = 0; j<pStdNames.length; j++) {
			pIndices.put(pStdNames[j], pPrior[j]);
		}
		//HashMap<String,String[]> synonyms = pAlphabet.getSynonyms();
		//get the transition matrices first, before calculating them over and over again
		double[][][] pMatrices = getMatrices(tree, pQmatrix);
		
		for(int i=0; i<cols.length; i++){
			String iCol = cols[i];
			for(int j=i+1; j<cols.length; j++){
				String jCol = cols[j];
				String[] colArray = {iCol, jCol};
				String[] leafData = toLeafData(colArray);
				BigDecimal likelihood = null;
				if(isNonstandardOverThreshold(leafData, pAlphabet, 0.5)) {
					likelihood = BigDecimal.ZERO;
				//} else if ((highlyConservedCols[i])||(highlyConservedCols[j])) {
				//	likelihood = pproduct(synonyms,pIndices,leafData);
				} else {
					likelihood = postOrderTraversal(tree, pQmatrix, pAlphabet, pPrior, leafData,pMatrices);
				}	
				probs.setPairingProbability(i, j, likelihood);				
			}
		}		
	
		return probs;
	}

	/*
	private boolean[] findConservedCols(String[] cols, Alphabet alphabet) {
		int length = cols.length;
		String[] standard = alphabet.getStandardNames();
		boolean[] conserved = new boolean[length];
		for (int i = 0; i<length; i++) {
			int ltrIdx = 0;
			while (ltrIdx<standard.length) {
				int count = 0;
				String[] tmp = cols[i].split("");
				//note empty string at start of array
				for (int j = 1; j<tmp.length; j++) {
					if (tmp[j].equals(standard[ltrIdx])) {
						count++;
					}
				}
				if (count==standard.length) {
					break;
				} else {
					ltrIdx++;
				}	
			}
			if (ltrIdx<standard.length) {
				conserved[i] = true;
			}
		}
		return conserved;
	}
	
	private BigDecimal pproduct(HashMap<String,String[]> synonyms, HashMap<String,Double> pIndices, String[] LeafData) {
		BigDecimal result = BigDecimal.ONE;
		for (int i = 0; i<LeafData.length; i++) {
			String[] syns = synonyms.get(LeafData[i]);
			double tmp = 0;
			for (String syn : syns){
				tmp += pIndices.get(syn);
			}
			result = result.multiply(BigDecimal.valueOf(tmp));
		}
		return result;
	}*/
	
	/**
	 * Helper method to calculate and return the transition probabilities matrix 
	 * Where Q is the rate matrix provided by the Knudsen-Hein model, Q = R*D*Rinverse
	 * @param R
	 * @param Rinverse
	 * @param D
	 * @return the matrix of transition probabilities
	 */
	private double[][] getTransitionProbsMatrix(Qmatrix qMatrix, double branchLengthDist){
		double D[][] = qMatrix.getDmatrix();
		double RI[][] = qMatrix.getRImatrix();
		double R[][] = qMatrix.getRmatrix();
		
		double[] eToDt = new double[D.length];
		for(int i=0; i<D.length; i++)
			eToDt[i] = Math.exp(D[i][i]*branchLengthDist);
		
		double[][] mtx = mult(multD(R,eToDt), RI);
		
		//R*e^Dt*R^-1			
		return mtx;
	}

	private double[][] multD(double[][] R, double[] ED) {
		double[][] result = new double[R.length][R.length];
		for (int i = 0; i<R.length; i++) {
			for (int j = 0; j<R.length; j++){
				result[i][j] = R[i][j]*ED[j];
			}
		}
		return result;
	}

	//Naive matrix multiplication -- would have used library if we needed more functions
	private double[][] mult(double[][] left, double[][] right){
		int leftRows = left.length,
			leftCols = left[0].length,
			rightCols = right[0].length;
		
		double[][] result = new double[leftRows][rightCols];
		
		  for(int i = 0; i < leftRows; i++) 
			    for(int j = 0; j < rightCols; j++) { 
			    	result[i][j] = 0;
			      for(int k = 0; k < leftCols; k++)
			        result[i][j] = result[i][j]+(left[i][k]*right[k][j]);
			    } 
		  
		return result;
	}

	/**
	 * Helper method to get columns of alignment
	 * @param alignment - String array representing a sequence alignment
	 * @return a string array where each string represents a column of the alignment
	 */
	private String[] convertToCols(String[] alignment)
	{
		String[] cols = new String[alignment[0].length()];

		for(int i=0; i<cols.length; i++)
			cols[i]="";
		for(int seqIdx=0; seqIdx<alignment.length; seqIdx++)
			for(int colIdx=0; colIdx<alignment[0].length(); colIdx++)
			{
				String currBase = alignment[seqIdx].substring(colIdx, colIdx+1);
				cols[colIdx] = cols[colIdx].concat(currBase);
			}
		return cols;
	}
	
	/**
	 * Helper method to get correct format of leafData
	 * @param cols - 1 or two strings representing 2 columns of an alignment
	 * @return an array of Strings (size = # seq in alignment) where each 
	 * string is data at a leaf
	 */
	private String[] toLeafData(String[] cols){
		String[] leafData = new String[cols[0].length()];
		for(int i=0; i<cols[0].length(); i++){
			//there are always exactly 1 or 2 columns
				leafData[i] = cols[0].substring(i, i+1).concat(cols[1].substring(i, i+1));		
		}
		return leafData;
	}
	
	private String[] toLeafData(String cols){
		String[] leafData = new String[cols.length()];
		for(int i=0; i<cols.length(); i++){
			leafData[i] = cols.substring(i, i+1);			
		}
		return leafData;
	}
	
	/**
	 * A method to determine if the number of non-standard letters in paired
	 * leafData is over some threshold.
	 * @param leafData - String[] representing two columns in the alignment
	 * @param alphabet - 
	 * @param threshold - 
	 * @return true if the number of nonstandard bases over threshold, and should 
	 * not allow pairs and false otherwise. 
	 */
	
	private boolean isNonstandardOverThreshold(String[] leafData, Alphabet alphabet, double threshold){
		int nonstandardCount = 0;
		double ratio =0;
		HashMap<String, String[]> synonyms = alphabet.getSynonyms();
		for(String s: leafData) {
			//nonstandard check
			if(synonyms.get(s).length > 1){
				nonstandardCount++;
			}
		}	
		ratio = (double) nonstandardCount/leafData.length;
		if(ratio > threshold)
			return true;

		else
			return false;	
	}
	
	
}



/*
//Naive matrix multiplication w/ BigDecimals 
private BigDecimal[][] mult(BigDecimal[][] left, BigDecimal[][] right){
	int leftRows = left.length,
		leftCols = left[0].length,
		rightCols = right[0].length;
	
	BigDecimal[][] result = new BigDecimal[leftRows][rightCols];
	
	  for(int i = 0; i < leftRows; i++) {
		    for(int j = 0; j < rightCols; j++) { 
		    	result[i][j] = BigDecimal.valueOf(0);
		      for(int k = 0; k < leftCols; k++) { 
		        result[i][j] = result[i][j].add(left[i][k].multiply(right[k][j]));
		      }
		    } 
		  }
	  
	return result;
}
*/
/*
//current node is not a leaf, so iterate through the node's children
//	for(String names1: standardLetters){
for(int ltrIdx1=0; ltrIdx1<standardLetters.length; ltrIdx1++){
	//BigDecimal tmp = BigDecimal.ONE;
	double tmp = 1.0;
	Enumeration children = nextNode.children();
	while(children.hasMoreElements()){
		DefaultMutableTreeNode nextChild = (DefaultMutableTreeNode) children.nextElement();
		EvolutionaryNodeData childData = (EvolutionaryNodeData) nextChild.getUserObject();
		int childIdx = childData.getIdx();
		//BigDecimal tmpp = BigDecimal.ZERO;
		double tmpp = 0.0;
		
		double distance = tree.getDistance(nextNode, nextChild);
		double[][] baseSubMatrix = getTransitionProbsMatrix(qMatrix, distance);
		
		for(int ltrIdx2=0; ltrIdx2<standardLetters.length; ltrIdx2++){
			//BigDecimal tmpTransProb = BigDecimal.valueOf(baseSubMatrix[standardLetters.indexOf(names1)][standardLetters.indexOf(names2)]);
			double tmpTransProb = baseSubMatrix[ltrIdx1][ltrIdx2];
			//tmpp = tmpp.add(L[traversal.indexOf(nextChild)][standardLetters.indexOf(names2)].multiply(tmpTransProb));
			tmpp = tmpp + (L[childIdx][ltrIdx2].doubleValue()* tmpTransProb);
		}
		//tmp = tmp.multiply(tmpp,MathContext.DECIMAL32);
		tmp = tmp * tmpp;
	}
	
	
	L[traversalIdx][ltrIdx1] = BigDecimal.valueOf(tmp);
}
} */