package uk.ac.ox.osscb.domain;

import javax.swing.tree.DefaultMutableTreeNode;

public class TreeNodeWithIdx {
	
	private DefaultMutableTreeNode treeNode;
	
	private Integer idx;
	
	public TreeNodeWithIdx(DefaultMutableTreeNode treeNode) {
		super();
		this.treeNode = treeNode;
	}

	public TreeNodeWithIdx(DefaultMutableTreeNode treeNode, int idx) {
		this(treeNode);
		setIdx(idx);
	}

	public DefaultMutableTreeNode getTreeNode() {
		return treeNode;
	}

	public int getIdx() {
		if(null == this.idx)
			throw new IllegalStateException("idx is not initialised");
		
		return this.idx.intValue();
	}

	public void setIdx(int idx) {
		if(null != this.idx){
			throw new IllegalStateException(String.format(
				"idx has already been initialised to: %d. Attempt to initialise %d",
				this.idx.intValue(), idx));
		}
		this.idx = idx;
	}
}
