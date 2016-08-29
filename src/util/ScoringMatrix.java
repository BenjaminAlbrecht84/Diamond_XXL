package util;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Vector;

public class ScoringMatrix {

	private String type;
	private int gapOpen, gapExtend;
	private Vector<String> alphabet;
	private int[][] matrix;

	public ScoringMatrix(String type, int gapOpen, int gapExtend) {
		this.type = type;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
//		parseMatrix(new File("./NCBI_ScoringMatrices/" + type));
		
		InputStream is  = this.getClass().getResourceAsStream("/"+type);
		parseMatrix(is);
		
	}

	public int[] cmpAlignmentScores(String s1, String s2) {
		boolean isGapOpen = false;
		int[] scores = new int[s1.length()];
		for (int i = 0; i < s1.length(); i++) {

			char a = s1.charAt(i);
			char b = s2.charAt(i);
			int score = -gapOpen - gapExtend;
			if (a == '-' || b == '-') {
				if (isGapOpen)
					score = -gapExtend;
				isGapOpen = true;
			} else {
				score = getScore(a, b);
				isGapOpen = false;
			}

			scores[i] = score;

		}
		return scores;
	}

	public int getScore(char a, char b) {
		int i = getIndex(a);
		int j = getIndex(b);
		return matrix[i][j];
	}

	private int getIndex(char c) {
		int i = (int) c;
		if (i >= 97 && i <= 122)
			return i - 65 - 32;
		if (i >= 65 && i <= 90)
			return i - 65;
		return 26;
	}

	private int parseMatrix(InputStream is) {

		try {

			alphabet = new Vector<String>();
			BufferedReader buf = new BufferedReader(new InputStreamReader(is));
			for (int i = 0; i < 2; i++) {
				String l = buf.readLine();
				Vector<String> entries = splitLine(l);
				alphabet.addAll(entries);
			}

			matrix = new int[27][27];
			int i = 0;
			String l;
			Vector<String> buffer = new Vector<String>();
			while ((l = buf.readLine()) != null) {
				if (!l.isEmpty()) {
					buffer.addAll(splitLine(l));
					if (buffer.size() == alphabet.size()) {
						int j = 0;
						int row = getIndex(alphabet.get(i).charAt(0));
						for (String e : buffer) {
							int val = Integer.parseInt(e);
							int col = getIndex(alphabet.get(j).charAt(0));
							matrix[row][col] = val;
							j++;
						}
						buffer.clear();
						i++;
					}
				}

			}
			buf.close();

			return 1;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return 0;
	}

	private Vector<String> splitLine(String l) {
		String[] entries = l.split("\\s+");
		Vector<String> splitResult = new Vector<String>();
		for (String e : entries) {
			if (!e.isEmpty() && !e.startsWith("/*") && !e.startsWith("*/"))
				splitResult.add(e.replaceAll(",", ""));
		}
		return splitResult;
	}

	public int getGapOpen() {
		return gapOpen;
	}

	public int getGapExtend() {
		return gapExtend;
	}

	public void print() {

		String out = type + " Matrix:\n";
		out = "\t";
		for (String c : alphabet)
			out = out.concat(c + "\t");
		out = out.concat("\n");

		for (int i = 0; i < alphabet.size(); i++) {
			out = out.concat(alphabet.get(i) + "\t");
			for (int j = 0; j < alphabet.size(); j++) {
				int row = getIndex(alphabet.get(i).charAt(0));
				int col = getIndex(alphabet.get(j).charAt(0));
				out = out.concat(matrix[row][col] + "\t");
			}
			out = out.concat("\n");
		}
		System.out.println(out);

	}

}
