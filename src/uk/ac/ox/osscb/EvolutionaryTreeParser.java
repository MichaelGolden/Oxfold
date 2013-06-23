package uk.ac.ox.osscb;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import javax.swing.tree.DefaultMutableTreeNode;


public class EvolutionaryTreeParser {
	
//
	/**
	 * Parses .newick file and puts data into an abstract evolutionary tree
	 * @author lmath
	 */

	/**
	 * parse method to produce tree
	 * @param path -- path to .newick file to be parsed
	 * @return an EvolutionaryTree parsed from .newick file
	 */
	public EvolutionaryTree parse(String path) {
		EvolutionaryTree tree = new EvolutionaryTree();
		try {
			Scanner sc = new Scanner(new File(path));
			if (sc.hasNextLine()) {
				String line = sc.nextLine().replaceAll("\n", "");			
				int strIdx = 1;
				tree = new EvolutionaryTree();
				EvolutionaryNodeData rootData = new EvolutionaryNodeData();
				tree.getRoot().setUserObject(rootData);
				//we'll iterate through the whole line char by char
				//line.length()-2 because we skip () containing the whole tree as well as 
				//terminating ;
				while(strIdx <= line.length()-2){
					boolean endSubtree = false;
					DefaultMutableTreeNode currentParent = tree.getRoot();
					DefaultMutableTreeNode node = null;
					// ; denotes the end of the tree in a newick file
					while(line.charAt(strIdx) != ';'){
						//we need to store the double value following any :
						if(line.charAt(strIdx)==':'){
							//double representation
							double currBranchLength = -1;
							//string representation, so we can change strIdx
							String currBranchLengthStr = null;
							
							//when we see the pattern ):someLength, this means that someLength belongs 
							//to the preceding subtree
							if(line.charAt(strIdx-1)==')')
								endSubtree = true;
							strIdx++;
							
							//get the double starting at current strIdx
							currBranchLengthStr = getDoubleStr(line, strIdx);
							strIdx += currBranchLengthStr.length();
							currBranchLength = new Double(currBranchLengthStr);
							
							//if this length belongs to the current node, just add the node to the 
							//subtree
							if(!endSubtree)
							{
						//		EvolutionaryNodeData data = new EvolutionaryNodeData(currBranchLength);
								node = new DefaultMutableTreeNode(new EvolutionaryNodeData(currBranchLength));
								currentParent.add(node);								
							}
							//otherwise add the branch length to the subtree's root
							else
							{
								((EvolutionaryNodeData) currentParent.getUserObject()).setBranchLength(currBranchLength);
							//	currentParent.setUserObject(new EvolutionaryNodeData(currBranchLength));
								DefaultMutableTreeNode parent = (DefaultMutableTreeNode) currentParent.getParent();
								currentParent = parent;
								endSubtree = false;						
							}
								
							strIdx++;
									
						}
						
						//If we see a (, create a new node to be parent to the following nodes
						//and assign the new node to be a child of the current parentNode
						else if(line.charAt(strIdx)=='(' && strIdx>0){
							DefaultMutableTreeNode parent = currentParent;
							currentParent = new DefaultMutableTreeNode(new EvolutionaryNodeData());
							parent.add(currentParent);
							strIdx++;
						}

					//we don't care what the nodes are named; just skip these characters.	
					else{
						strIdx++;
					}											
						
					}
				}
				sc.close();
			}
			
		} catch (FileNotFoundException e) {
			throw new IllegalArgumentException(String.format("File %s not found.",path));
		}
		return tree;
	}
	
	/**
	 * Helper method to get a string representation of a double.
	 * @param line -- whole line representing Newick structure
	 * @param idx -- where the double starts
	 * @return a string representation of the double starting at idx.
	 */
	private String getDoubleStr(String line, int idx){
		StringBuilder sb = new StringBuilder();
		while(line.substring(idx, idx+1).matches("[0-9]*[.]?[0-9]*")){
			sb.append(line.charAt(idx));
			idx++;
		}
	 
		return sb.toString();
	}
	


}



