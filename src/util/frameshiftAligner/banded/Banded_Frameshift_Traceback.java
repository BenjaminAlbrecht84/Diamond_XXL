package util.frameshiftAligner.banded;

import java.util.HashMap;
import java.util.Vector;

import util.ScoringMatrix;
import util.frameshiftAligner.banded.Banded_Frameshift_Alignment.AliMode;

public class Banded_Frameshift_Traceback {

	private ScoringMatrix scoringMatrix;
	private int gapOpen, gapExtend, delta;

	public Banded_Frameshift_Traceback(ScoringMatrix scoringMatrix, int gapOpen, int gapExtend, int delta) {
		this.scoringMatrix = scoringMatrix;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.delta = delta;
	}

	public Object[] run(int[][] D1, int[][] Q1, int[][] P1, int[][] D2, int[][] Q2, int[][] P2, int[][] D3, int[][] Q3, int[][] P3, String s1,
			String s11, String s12, String s13, String s2, int[] startingCell, Vector<int[]> bounderies, AliMode mode) {

		HashMap<Integer, int[][][]> frameToMatrix = new HashMap<Integer, int[][][]>();
		int[][][] frame1Matrices = { D1, Q1, P1 };
		frameToMatrix.put(1, frame1Matrices);
		int[][][] frame2Matrices = { D2, Q2, P2 };
		frameToMatrix.put(2, frame2Matrices);
		int[][][] frame3Matrices = { D3, Q3, P3 };
		frameToMatrix.put(3, frame3Matrices);

		HashMap<Integer, String> frameToSequence = new HashMap<Integer, String>();
		frameToSequence.put(1, s11);
		frameToSequence.put(2, s12);
		frameToSequence.put(3, s13);

		StringBuffer ali1 = new StringBuffer();
		StringBuffer ali2 = new StringBuffer();
		StringBuffer frames = new StringBuffer();
		StringBuffer[] alignment = { ali1, ali2, frames };

		// running traceback
		int[] cell = startingCell;
		while (true) {

			cell = backtrace(frameToMatrix, frameToSequence, s2, cell, alignment, bounderies, mode);

			int i = cell[0];
			int j = cell[1];

			if (bounderies != null && (j == bounderies.get(i)[0] || j == bounderies.get(i)[1]))
				return null;
			else if (i == 0 && j == 0)
				break;
			else if (j == 0 && mode == AliMode.SEMI_GLOBAL)
				break;

		}

		// reversing sequences
		if (mode == AliMode.FREESHIFT_LEFT) {
			alignment[0] = alignment[0].reverse();
			alignment[1] = alignment[1].reverse();
			alignment[2] = alignment[2].reverse();
		}

		int[] alignedInfo = mode == AliMode.FREESHIFT_LEFT ? startingCell : mode == AliMode.FREESHIFT_RIGHT ? startingCell : cell;
		int aliScore = startingCell[3];
		Object[] res = { alignment[0].toString(), alignment[1].toString(), alignment[2].toString(), aliScore, alignedInfo, s11.length() };
		return res;

	}

	private int[] backtrace(HashMap<Integer, int[][][]> frameToMatrix, HashMap<Integer, String> frameToSequence, String s2, int[] cell,
			StringBuffer[] alignment, Vector<int[]> bounderies, AliMode mode) {

		int i = cell[0];
		int j = cell[1];
		int frame = cell[2];

		int[] frames = { 1, 2, 3 };
		switch (frame) {
		case 2:
			frames[0] = 2;
			frames[1] = 1;
			break;
		case 3:
			frames[0] = 3;
			frames[2] = 1;
			break;
		}

		int[][] D1 = frameToMatrix.get(frames[0])[0];
		int[][] Q1 = frameToMatrix.get(frames[0])[1];
		int[][] P1 = frameToMatrix.get(frames[0])[2];
		String s1 = frameToSequence.get(frames[0]);
		int[] prevCell = null;

		int penalty = 0;
		for (int f : frames) {

			// searching previous cell
			int[][] D2 = frameToMatrix.get(f)[0];
			int[][] Q2 = frameToMatrix.get(f)[1];
			int[][] P2 = frameToMatrix.get(f)[2];
			prevCell = findPreviousCell(D1, Q1, P1, D2, Q2, P2, i, j, alignment, s1, s2, f, frames[0], penalty, frameToSequence);
			if (prevCell != null)
				break;

			// setting penalty for changing frames
			penalty = delta;

		}

		return prevCell;

	}

	private int[] findPreviousCell(int[][] D1, int[][] Q1, int[][] P1, int[][] D2, int[][] Q2, int[][] P2, int i, int j, StringBuffer[] alignment,
			String s1, String s2, int frame, int curFrame, int penalty, HashMap<Integer, String> frameToSequence) {

		if (i != 0 || j != 0) {

			// upper border reached
			if (i == 0) {
				StringBuffer subSeq = new StringBuffer("");
				for (int k = j - 1; k >= 0; k--)
					subSeq.insert(0, s2.charAt(k));
				alignment[0] = alignment[0].insert(0, gapString(j));
				alignment[1] = alignment[1].insert(0, subSeq);
				alignment[2] = alignment[2].insert(0, frameString(j, curFrame));
				int[] res = { i, 0, frame };
				return res;
			}

			// left border reached
			else if (j == 0) {
				StringBuffer subSeq = new StringBuffer("");
				for (int k = i - 1; k >= 0; k--)
					subSeq.insert(0, s1.charAt(k));
				alignment[0] = alignment[0].insert(0, subSeq);
				alignment[1] = alignment[1].insert(0, gapString(i));
				alignment[2] = alignment[2].insert(0, frameString(i, curFrame));
				int[] res = { 0, j, frame };
				return res;
			}

			else {

				// checking diagonal
				if (i != 0 && j != 0 && D2[i - 1][j - 1] + scoringMatrix.getScore(s1.charAt(i - 1), s2.charAt(j - 1)) - penalty == D1[i][j]) {
					alignment[0] = alignment[0].insert(0, s1.charAt(i - 1));
					alignment[1] = alignment[1].insert(0, s2.charAt(j - 1));
					alignment[2] = alignment[2].insert(0, curFrame);
					int[] res = { i - 1, j - 1, frame };
					return res;
				}

				else {

					// checking left direction
					Integer l = null;
					int score = D1[i][j];
					Vector<Integer> frames = new Vector<Integer>();
					for (int k = j - 1; k >= 0; k--) {

						if (D2[i][k] == Integer.MIN_VALUE)
							break;
						if (D2[i][k] - w(j - k) - penalty == D1[i][j]) {
							frames.add(0, frame);
							l = k;
							break;
						}

						score += gapExtend;
						if (Q1[i][k] != score && Q2[i][k] - penalty != score)
							break;
						if (Q1[i][k] == score)
							frames.add(0, curFrame);
						else
							frames.add(0, frame);

					}
					if (l != null) {
						StringBuffer subSeq = new StringBuffer("");
						for (int k = j - 1; k >= l; k--)
							subSeq.insert(0, s2.charAt(k));
						alignment[0] = alignment[0].insert(0, gapString(j - l));
						alignment[1] = alignment[1].insert(0, subSeq);
						// alignment[2] = alignment[2].insert(0, frameString(j - l, curFrame));
						alignment[2] = alignment[2].insert(0, frameString(frames));
						int[] res = { i, l, frame };
						return res;
					}

					else {

						// checking upper direction
						Integer t = null;
						score = D1[i][j];
						frames = new Vector<Integer>();
						for (int k = i - 1; k >= 0; k--) {
							if (D2[k][j] == Integer.MIN_VALUE)
								break;
							if (D2[k][j] - w(i - k) - penalty == D1[i][j]) {
								frames.add(0, frame);
								t = k;
								break;
							}
							score += gapExtend;
							if (P1[k][j] != score && P2[k][j] - penalty != score)
								break;
							if (P1[k][j] == score)
								frames.add(0, curFrame);
							else
								frames.add(0, frame);

						}
						if (t != null) {
							StringBuffer subSeq = new StringBuffer("");
							// for (int k = i - 1; k >= t; k--)
							// subSeq.insert(0, s1.charAt(k));
							int pos = frames.size() - 1;
							for (int k = i - 1; k >= t; k--) {
								String s = frameToSequence.get(frames.get(pos--));
								subSeq.insert(0, s.charAt(k));
							}
							alignment[0] = alignment[0].insert(0, subSeq);
							alignment[1] = alignment[1].insert(0, gapString(i - t));
							// alignment[2] = alignment[2].insert(0, frameString(i - t, curFrame));
							alignment[2] = alignment[2].insert(0, frameString(frames));
							int[] res = { t, j, frame };
							return res;
						}

					}
				}
			}

		}

		return null;

	}

	private Object gapString(Integer l) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < l; i++)
			buf.append("-");
		return buf.toString();
	}

	private Object frameString(Integer l, int frame) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < l; i++)
			buf.append(frame);
		return buf.toString();
	}

	private Object frameString(Vector<Integer> frames) {
		StringBuffer buf = new StringBuffer();
		for (int f : frames)
			buf.append(f);
		return buf.toString();
	}

	private int w(int k) {
		return gapExtend * k + gapOpen;
	}

}
