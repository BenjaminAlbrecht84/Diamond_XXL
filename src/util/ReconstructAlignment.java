package util;

import java.util.BitSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReconstructAlignment {

	private ScoringMatrix matrix;

	public ReconstructAlignment(ScoringMatrix matrix) {
		this.matrix = matrix;
	}

	public Object[] run(String seq, String cigar, String md) {

		// System.out.println("\n>" + seq + " " + cigar + " " + md + "\n");

		// calculate query string
		Object[] queryResult = reconstructQuerySeq(seq, cigar);
		String query = (String) queryResult[0];
		Vector<Integer[]> insertions = (Vector<Integer[]>) queryResult[1];
		for (int i = insertions.size() - 1; i >= 0; i--) {
			Integer[] pair = insertions.get(i);
			seq = new StringBuilder(seq).replace(pair[0], pair[0] + pair[1], "").toString();
		}

		// calculate ref string
		String ref = reconstructRefSeq(seq, md, (Vector<Integer[]>) queryResult[1], query);

		String[] alignments = { query, ref };
		return run(alignments);
	}

	public Object[] run(String[] alignments) {

		String query = alignments[0];
		String ref = alignments[1];

		int[] scores = matrix.cmpAlignmentScores(query, ref);
		int score = 0;
		for (int s : scores)
			score += s;

		BitSet query_insertions = new BitSet(scores.length);
		for (int i = 0; i < ref.length(); i++) {
			if (ref.charAt(i) == '-')
				query_insertions.set(i);
		}

		BitSet query_deletions = new BitSet(scores.length);
		for (int i = 0; i < query.length(); i++) {
			if (query.charAt(i) == '-')
				query_deletions.set(i);
		}

		Object[] result = { scores, query_insertions, query_deletions, score, query, ref };

		return result;
	}

	private Object[] reconstructQuerySeq(String seq, String cigar) {

		Matcher matcher = Pattern.compile("[0-9]+[MID]+").matcher(cigar);
		Vector<String> pairs = new Vector<String>();
		while (matcher.find())
			pairs.add(matcher.group());

		StringBuffer query = new StringBuffer("");
		int i = 0, d = 0;
		Vector<Integer[]> insertions = new Vector<Integer[]>();
		for (String p : pairs) {
			int num = Integer.parseInt(p.substring(0, p.length() - 1));
			char type = p.charAt(p.length() - 1);
			if (type == 'M' || type == 'I') {
				query = query.append(seq.substring(i, i + num));
				if (type == 'I') {
					// System.out.println(query.length() - num - d);
					Integer[] ins = { query.length() - num - d, num };
					insertions.add(ins);
				}
				i += num;
			} else if (type == 'D') {
				String gaps = "";
				for (int k = 0; k < num; k++)
					gaps = gaps.concat("-");
				d += num;
				query = query.append(gaps);
			} else
				System.err.print("ERROR: Parsing CIGAR-String - unknown type " + type);
		}

		Object[] result = { query.toString(), insertions };
		return result;
	}

	private String reconstructRefSeq(String seq, String md, Vector<Integer[]> insertions, String query) {

		Vector<MD_Item> items = new Vector<MD_Item>();
		MD_Item curr_item = new MD_Item(md.charAt(0));
		for (int i = 1; i < md.length(); i++) {
			char c = md.charAt(i);
			if (typeOf(curr_item) != typeOf(c)) {
				items.add(curr_item);
				curr_item = new MD_Item(md.charAt(i));
			} else
				curr_item.addCharacter(c);
		}
		items.add(curr_item);

		StringBuffer refBuf = new StringBuffer("");
		int pos = 0;
		Vector<Integer> deletionSites = new Vector<Integer>();
		for (MD_Item item : items) {
			if (item.getType() == 'M') {
				int num = Integer.parseInt(item.getSeq());
				refBuf = refBuf.append(seq.substring(pos, pos + num));
				pos += num;
			} else if (item.getType() == 'X') {
				refBuf = refBuf.append(item.getSeq());
				pos += item.getSeq().length();
			} else if (item.getType() == 'D') {
				for (int i = refBuf.length(); i < refBuf.length() + item.getSeq().length(); i++) {
					deletionSites.add(i);
				}
				refBuf = refBuf.append(item.getSeq());
			} else
				System.err.print("ERROR: Parsing MD-String - unknown type " + item.getType());
		}

		BitSet deletionSet = new BitSet(refBuf.length());
		for (int i = 0; i < refBuf.length(); i++)
			deletionSet.set(i, deletionSites.contains(i));

		int inserted = 0;
		for (int i = 0; i < insertions.size(); i++) {

			int p = insertions.get(i)[0];
			int num = insertions.get(i)[1];

			for (int k = 0; k < p - inserted; k++) {
				if (deletionSet.get(k))
					p++;
			}
			while (query.charAt(p) == '-')
				p++;

			String gaps = "";
			for (int k = 0; k < num; k++)
				gaps = gaps.concat("-");
			refBuf = refBuf.insert(p, gaps);
			inserted += num;
		}

		return refBuf.toString();

	}

	private char typeOf(MD_Item item) {
		if (!item.getSeq().isEmpty())
			return typeOf(item.getSeq().charAt(0));
		return 'C';
	}

	private char typeOf(char c) {
		int i = (int) c;
		if (i >= 48 && i <= 57)
			return 'N';
		if (c == '^')
			return 'D';
		return 'C';
	}

	private char getItemType(char c) {
		int i = (int) c;
		if (i >= 48 && i <= 57)
			return 'M';
		else if (c == '^')
			return 'D';
		return 'X';
	}

	public class MD_Item {

		private char type;
		private StringBuffer seq = new StringBuffer("");

		public MD_Item(char c) {
			this.type = getItemType(c);
			if (type != 'D')
				addCharacter(c);
		}

		public void addCharacter(char c) {
			seq = seq.append(c);
		}

		public char getType() {
			return type;
		}

		public String getSeq() {
			return seq.toString();
		}

	}

}
