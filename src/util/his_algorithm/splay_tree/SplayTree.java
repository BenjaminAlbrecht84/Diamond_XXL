package util.his_algorithm.splay_tree;

public class SplayTree {

	private SplayNode root;

	public SplayNode getMin() {
		SplayNode v = root;
		if (v == null)
			return null;
		while (v.getLeftChild() != null)
			v = v.getLeftChild();
		splay(v);
		return v;
	}

	public SplayNode prevNode(int id) {
		SplayNode r = null;
		if (root != null)
			r = prevNodeRec(root, id, r);
		splay(r);
		return r;
	}

	private SplayNode prevNodeRec(SplayNode v, int id, SplayNode r) {
		if (v != null) {
			if (r == null && v.getId() < id)
				r = v;
			else if (r != null && r.getId() < v.getId() && v.getId() < id)
				r = v;
			if (v.getId() >= id)
				return prevNodeRec(v.getLeftChild(), id, r);
			return prevNodeRec(v.getRightChild(), id, r);
		}
		return r;
	}

	public SplayNode nextNode(int id) {
		SplayNode r = null;
		if (root != null)
			r = nextNodeRec(root, id, r);
		splay(r);
		return r;
	}

	private SplayNode nextNodeRec(SplayNode v, int id, SplayNode r) {
		if (v != null) {
			if (r == null && v.getId() > id)
				r = v;
			else if (r != null && r.getId() > v.getId() && v.getId() > id)
				r = v;
			if (v.getId() >= id)
				return nextNodeRec(v.getRightChild(), id, r);
			return nextNodeRec(v.getLeftChild(), id, r);
		}
		return r;
	}

	public SplayNode insert(int id, double weight) {
		SplayNode v = new SplayNode(id, weight);
		if (root == null)
			root = v;
		else {
			insertNodeRec(root, v);
			splay(v);
		}
		return v;
	}

	private void insertNodeRec(SplayNode p, SplayNode v) {
		if (p.isLeaf()) {
			p.addChild(v);
		} else {
			SplayNode c = p.getChild(v.getId());
			if (c == null) {
				p.addChild(v);
			} else
				insertNodeRec(c, v);
		}
	}

	public void delete(SplayNode v) {
		SplayNode p = v.getParent();
		if (v.getNumOfChildren() == 0) {
			if (p == null)
				root = null;
			else
				p.removeChild(v);
		} else if (v.getNumOfChildren() == 1) {
			supressNode(v);
		} else {
			SplayNode w = v.getRightChild();
			while (w.getLeftChild() != null)
				w = w.getLeftChild();
			v.copyNode(w);
			if (w.getNumOfChildren() == 1)
				supressNode(w);
			else
				w.getParent().removeChild(w);
		}
		if (p != null)
			splay(p);
	}

	private void supressNode(SplayNode v) {
		SplayNode p = v.getParent();
		SplayNode c = v.getLeftChild() != null ? v.getLeftChild() : v.getRightChild();
		if (p == null) {
			setRoot(c);
		} else {
			boolean vLeft = p.getLeftChild() == null ? false : p.getLeftChild().equals(v);
			if (vLeft)
				p.setLeftChild(c);
			else
				p.setRightChild(c);
		}
	}

	// different splay-cases mentioned in wikipedia
	// -> see https://en.wikipedia.org/wiki/Splay_tree
	private void splay(SplayNode x) {

		// print();

		if (x != null && !x.equals(root)) {
			SplayNode p = x.getParent();
			boolean xLeft = p.getLeftChild() == null ? false : p.getLeftChild().equals(x);
			if (p.equals(root)) {
				p.removeChild(x);
				if (xLeft) {
					// System.out.println("Case 0.L");
					SplayNode r = x.getRightChild();
					setRoot(x);
					x.setRightChild(p);
					if (r != null)
						p.setLeftChild(r);
				} else {
					// System.out.println("Case 0.R");
					SplayNode l = x.getLeftChild();
					setRoot(x);
					x.setLeftChild(p);
					if (l != null)
						p.setRightChild(l);
				}
			} else {
				SplayNode g = p.getParent();
				boolean pLeft = g.getLeftChild() == null ? false : g.getLeftChild().equals(p);
				if (pLeft && xLeft) {
					// System.out.println("Case 1.L");
					SplayNode r = x.getRightChild();
					p.removeChild(x);
					flipLeftRec(p);
					setRoot(x);
					x.setRightChild(p);
					if (r != null)
						p.setLeftChild(r);
				} else if (!pLeft && !xLeft) {
					// System.out.println("Case 1.R");
					SplayNode l = x.getLeftChild();
					p.removeChild(x);
					flipRightRec(p);
					setRoot(x);
					x.setLeftChild(p);
					if (l != null)
						p.setRightChild(l);
				} else if (pLeft && !xLeft) {
					// System.out.println("Case 2.L");
					SplayNode r = x.getRightChild();
					g.removeChild(p);
					if (r != null)
						g.setLeftChild(r);
					flipLeftRec(g);
					SplayNode l = x.getLeftChild();
					p.removeChild(x);
					setRoot(x);
					x.setLeftChild(p);
					x.setRightChild(g);
					if (l != null)
						p.setRightChild(l);
				} else {
					// System.out.println("Case 2.R");
					SplayNode l = x.getLeftChild();
					g.removeChild(p);
					if (l != null)
						g.setRightChild(l);
					flipRightRec(g);
					SplayNode r = x.getRightChild();
					p.removeChild(x);
					setRoot(x);
					x.setLeftChild(g);
					x.setRightChild(p);
					if (r != null)
						p.setLeftChild(r);
				}
			}

		}

		// System.out.println();
		// print();
		// System.out.println("**********************");

	}

	private void flipLeftRec(SplayNode v) {
		SplayNode p = v.getParent();
		if (p != null) {
			flipLeftRec(p);
			SplayNode r = v.getRightChild();
			p.removeChild(v);
			if (r != null)
				p.setLeftChild(r);
			v.setRightChild(p);
		}
	}

	private void flipRightRec(SplayNode v) {
		SplayNode p = v.getParent();
		if (p != null) {
			flipRightRec(p);
			SplayNode l = v.getLeftChild();
			p.removeChild(v);
			if (l != null)
				p.setRightChild(l);
			v.setLeftChild(p);
		}
	}

	private void setRoot(SplayNode x) {
		x.setParent(null);
		root = x;
	}

	public SplayNode getRoot() {
		return root;
	}

	// ****************************************************

	public void print() {
		if (root != null)
			root.print();
	}

}
