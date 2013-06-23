package uk.ac.ox.osscb;
import java.util.Enumeration;

import javax.swing.tree.DefaultMutableTreeNode;

/**
 * Class to represent an evolutionary tree with decimal value branch lengths.
 * @author lmath@cs.ubc.ca
 *
 */

public class EvolutionaryTree {

	DefaultMutableTreeNode root;
	
	/**
	 * Constructor; children are added later by EvolutionaryTreeParser
	 */
	public EvolutionaryTree(){
		root = new DefaultMutableTreeNode();
	}
	
	/**
	 * Method for the total distance between two nodes.
	 * @param ancestor -- a node in the tree, an ancestor to descendant.
	 * @param descendant -- a node in the tree, a descendant of ancestor.
	 * @return distance between startNode and targetNode
	 */
	public double getDistance(DefaultMutableTreeNode ancestor, DefaultMutableTreeNode descendant){
		if(descendant.isNodeDescendant(ancestor))
			throw new IllegalArgumentException("One node must be a descendant of the other");
		
		DefaultMutableTreeNode parent = descendant;
		double totalDistance = 0;
		
		while(!parent.equals(ancestor)){
			EvolutionaryNodeData data = (EvolutionaryNodeData) parent.getUserObject();
			totalDistance+= data.getBranchLength();
			parent = (DefaultMutableTreeNode) parent.getParent();
		}
		
		return totalDistance;
	}
	
	/**
	 * Getter for root of tree
	 * @return root - the root of the tree
	 */
	public DefaultMutableTreeNode getRoot(){
		return root;
	}
	
	/**
	 * Method to get the total number of nodes in a tree rooted at node
	 * @param node - root
	 * @return total number of nodes in the tree
	 */
	public int nodesInTreeCount(DefaultMutableTreeNode node)
	{		
		if(node.isLeaf())
			return 1;
		
		else{
			int count = 0;
			Enumeration e = node.children();
			while (e.hasMoreElements())
				count += nodesInTreeCount((DefaultMutableTreeNode) e.nextElement());
		
			return 1 + count;
		}
	
	}
	
}
