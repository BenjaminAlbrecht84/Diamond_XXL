package util.his_algorithm.splay_tree;

public class SplayNode {

	private SplayNode parent;
	private SplayNode[] children = new SplayNode[2];
	private int id;
	private double weight;
	
	private Object[] info;

	public SplayNode(SplayNode parent, int id, double weight) {
		this.parent = parent;
		this.id = id;
		this.weight = weight;
	}

	public SplayNode(int id, double weight) {
		this.id = id;
		this.weight = weight;
	}

	public void addChild(SplayNode c) {
		if (c.getId() <= id)
			children[0] = c;
		else
			children[1] = c;
		c.setParent(this);
	}
	
	public void setLeftChild(SplayNode c) {
		children[0] = c;
		c.setParent(this);
	}

	public void setRightChild(SplayNode c) {
		children[1] = c;
		c.setParent(this);
	}

	public SplayNode getChild(int val) {
		if (val <= id)
			return children[0];
		return children[1];
	}

	public SplayNode getLeftChild() {
		return children[0];
	}

	public SplayNode getRightChild() {
		return children[1];
	}

	public boolean isLeaf() {
		if (children[0] == null && children[1] == null)
			return true;
		return false;
	}

	public SplayNode getParent() {
		return parent;
	}

	public void setParent(SplayNode parent) {
		this.parent = parent;
	}

	public int getId() {
		return id;
	}

	public void copyNode(SplayNode v) {
		this.id = v.getId();
		this.weight = v.getWeight();
	}

	public double getWeight() {
		return weight;
	}

	public int getNumOfChildren() {
		int counter = 0;
		for (SplayNode c : children) {
			if (c != null)
				counter++;
		}
		return counter;
	}

	public void removeChild(SplayNode w) {
		for (int i = 0; i < 2; i++) {
			if (children[i] == w)
				children[i] = null;
		}
	}	
	
	public Object[] getInfo() {
		return info;
	}

	public void setInfo(Object[] info) {
		this.info = info;
	}

	// *************************************************************

	public void print() {
		print("", true);
	}

	private void print(String prefix, boolean isTail) {
		if (children[1] != null)
			children[1].print(prefix + (isTail ? "    " : "    "), false);
		System.out.println(prefix + (isTail ? " ---- " : "+--- ") + id);
		if (children[0] != null)
			children[0].print(prefix + (isTail ? "    " : "    "), false);
	}

}
