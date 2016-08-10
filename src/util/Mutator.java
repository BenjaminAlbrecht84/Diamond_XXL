package util;

import java.util.Vector;

public class Mutator {

	private Character[] sigma;

	public Mutator(Character[] sigma) {
		this.sigma = sigma;
	}

	public String insert_deletions(String seq, double error_rate) {
		Vector<Integer> positions = new Vector<Integer>();
		for (int i = seq.length() - 1; i >= 0; i--) {
			if (Math.random() < error_rate)
				positions.add(i);
		}
		for (int pos : positions) {
				seq = deletePos(seq, pos);
		}
		return seq;
	}

	private String deletePos(String seq, int pos) {
		return new StringBuilder(seq).replace(pos, pos + 1, "").toString();
	}

	public String insert_insertions(String seq, double error_rate) {
		Vector<Integer> positions = new Vector<Integer>();
		for (int i = seq.length() - 1; i >= 0; i--) {
			if (Math.random() < error_rate)
				positions.add(i);
		}
		String mut_seq = seq;
		for (int p : positions)
			mut_seq = insertChar(mut_seq, p);
		return mut_seq;
	}

	private String insertChar(String seq, int pos) {
		return new StringBuilder(seq).insert(pos, getRandomChar()).toString();
	}

	public String insert_substitutions(String seq, double error_rate) {
		int length = seq.length();
		for (int i = 0; i < length; i++) {
			if (Math.random() < error_rate) 
				seq = substituteChar(seq,i);
		}
		return seq;
	}

	private String substituteChar(String seq, int pos) {
		char c = getRandomChar();
		while (c == seq.charAt(pos))
			c = getRandomChar();
		return new StringBuilder(seq).replace(pos, pos + 1, c + "").toString();
	}

	private char getRandomChar() {
		double l = sigma.length;
		int j = (int) Math.floor(Math.random() * l);
		j = (j == sigma.length) ? sigma.length - 1 : j;
		return sigma[j];
	}

}
