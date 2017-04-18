package util.frameshiftAligner.banded.Zheng;

import java.util.Vector;

import util.CodonTranslator_Array;
import util.ScoringMatrix;
import util.frameshiftAligner.xDrop.Zheng_XDrop_Frameshift.AliMode;

public class Zheng_Banded_Frameshift {

	public enum AliMode {
		GLOBAL, FREESHIFT_LEFT, FREESHIFT_RIGHT, SEMI_GLOBAL
	};

	public enum BORDER_REACHED {
		LEFT, RIGHT
	};

	private BORDER_REACHED borderReached;

	private CodonTranslator_Array codonTranslator = new CodonTranslator_Array();
	final static int MIN_VALUE = Integer.MIN_VALUE;
	private String protein, dna;
	private Vector<int[]> diagonal;

	private boolean verbose = false;

	private ScoringMatrix scoringMatrix;
	private int gop, gep, F;

	public Zheng_Banded_Frameshift(ScoringMatrix scoringMatrix, int delta) {
		this.scoringMatrix = scoringMatrix;
		this.gop = scoringMatrix.getGapOpen() + scoringMatrix.getGapExtend();
		this.gep = scoringMatrix.getGapExtend();
		this.F = delta;
	}

	public Object[] run(String s1, String s2, AliMode mode, double bandLeft, double bandRight) {

		// System.out.println(s1 + " " + s2);

		this.dna = s1.replace('T', 'U');
		this.protein = s2;

		if (mode == AliMode.FREESHIFT_LEFT) {
			protein = new StringBuilder(protein).reverse().toString();
			dna = new StringBuilder(dna).reverse().toString();
		}

		int n = protein.length() + 1;
		int m = dna.length() + 4;

		// initializing matrices
		int[][] X = new int[n][m];
		int[][] Y = new int[2][m];
		initMatrices(n, m, X, Y, mode);

		// computing diagonal
		Vector<int[]> tmpDiagonal = new Vector<int[]>();
		if (diagonal == null) {
			diagonal = new Vector<int[]>();
			tmpDiagonal = cmpDiagonal(X.length, X[0].length - 6);
			for (int[] p : tmpDiagonal) {
				int[] point = { p[0] + 1, p[1] + 4 };
				diagonal.add(point);
			}
		}

		// initializing diagonal pointer
		int pointer = 0;
		while (diagonal.get(pointer)[0] == 0)
			pointer++;

		// initializing diagonal bounds
		Vector<int[]> bounds = new Vector<int[]>();
		int[] initBound = { -1, m - 2 + 1 };
		bounds.add(initBound);

		// computing bandwidth
		int dLeft = computeBandwith(n, m - 2, bandLeft);
		int dRight = computeBandwith(n, m - 2, bandRight);

		// filling Matrix
		fillingMatrix(n, m, pointer, X, Y, mode, bounds, dLeft, dRight);
		if (verbose)
			printMatrix(X, protein, dna);

		// starting traceback
		Object[] traceback_result = traceback(X, mode, bounds);

		// return traceback_result;

		// returning result
		if (traceback_result == null) {
			if (borderReached == BORDER_REACHED.LEFT)
				bandLeft = 2. * bandLeft <= 1. ? 2. * bandLeft : 1.;
			if (borderReached == BORDER_REACHED.RIGHT)
				bandRight = 2. * bandRight <= 1. ? 2. * bandRight : 1.;
			s2 = mode == AliMode.FREESHIFT_LEFT ? new StringBuffer(s2).reverse().toString() : s2;
			return run(s1, s2, mode, bandLeft, bandRight);
		} else {
			diagonal = null;

			// printResult(traceback_result);

			return traceback_result;
		}

	}

	private void fillingMatrix(int n, int m, int pointer, int[][] X, int[][] Y, AliMode mode, Vector<int[]> bounds, int dLeft, int dRight) {

		for (int i = 1; i < n; i++) {

			// initializing band bounds
			int[] leftCell = diagonal.get(pointer);
			while (pointer < diagonal.size() && diagonal.get(pointer)[0] == i)
				pointer++;
			int[] rightCell = diagonal.get(pointer - 1);

			// adapting & storing band bounds
			int l = leftCell[1] - dLeft;
			int r = mode == AliMode.GLOBAL ? rightCell[1] + dRight : m - 2;
			int[] bound = { l, r };
			bounds.add(bound);

			// initializing Y
			Y[i % 2][0] = MIN_VALUE;
			for (int j = 1; j < 4; j++)
				Y[i % 2][j] = -gop - ((i - 1) * gep);
			for (int j = 4; j < m - 2; j++)
				Y[i % 2][j] = MIN_VALUE;

			// initializing D
			int[] D = { Y[i % 2][2] - gop, Y[i % 2][1] - gop, Y[i % 2][3] - gop };
			if (l >= 5)
				D[2] = MIN_VALUE;
			if (l >= 6)
				D[0] = MIN_VALUE;
			if (l >= 7)
				D[1] = MIN_VALUE;

			r = r > m - 2 ? m - 2 : r;
			l = l > 4 ? l : 4;
			for (int j = l; j < r; j++) {

				int yMax = Math.max(sub(X[i - 1][j], gop), sub(Y[(i - 1) % 2][j], gep));
				Y[i % 2][j] = yMax;

				D[j % 3] = Math.max(sub(X[i][j - 3], gop), sub(D[j % 3], gep));
				int dMax = Math.max(yMax, D[j % 3]);

				int a = X[i - 1][j - 3];
				if (j >= 6)
					a = Math.max(a, sub(X[i - 1][j - 2], F));
				if (j >= 8)
					a = Math.max(a, sub(X[i - 1][j - 4], F));

				int b = add(a, matchScore(i, j, false, mode));
				int xMax = Math.max(b, dMax);
				X[i][j] = xMax;

			}
		}
	}

	private int computeBandwith(int n1, int n2, double band) {
		int d = (int) Math.ceil((double) Math.max(n1, n2) * band);
		int a = n1 < n2 ? n1 : n2;
		int b = n2 >= n1 ? n2 : n1;
		double frac = (1. - ((double) a / (double) b));
		int dMin = (int) Math.floor((double) Math.max(n1, n2) * frac);
		d += dMin;
		return d;
	}

	private int[] cmpStartingCell(int[][] X, StringBuffer[] alignment, AliMode mode) {

		int max = Integer.MIN_VALUE;
		int[] startingCell = new int[3];

		// searching for max in last row
		for (int i = 3; i < X[0].length - 1; i++) {
			if (mode == AliMode.SEMI_GLOBAL && X[0].length - i <= 3)
				break;
			int v = X[X.length - 1][X[0].length - i];
			if (v >= max) {
				max = v;
				startingCell[0] = X.length - 1;
				startingCell[1] = X[0].length - i;
				startingCell[2] = v;
			}
			if (mode == AliMode.GLOBAL && i == 5)
				break;
		}

		if (mode == AliMode.SEMI_GLOBAL) {
			// searching for max in last column
			for (int j = 0; j < X.length; j++) {
				if (j == 0)
					break;
				int v = X[j][X[0].length - 3];
				if (v >= max) {
					max = v;
					startingCell[0] = j;
					startingCell[1] = X[0].length - 3;
					startingCell[2] = v;
				}
			}
		}

		return startingCell;

	}

	private Object[] traceback(int[][] X, AliMode mode, Vector<int[]> bounds) {

		StringBuffer ali1 = new StringBuffer();
		StringBuffer ali2 = new StringBuffer();
		StringBuffer frames = new StringBuffer();
		StringBuffer[] alignment = { ali1, ali2, frames };

		int[] startingCell = cmpStartingCell(X, alignment, mode);
		int[] cell = startingCell;

		borderReached = null;
		int lastCol = 0, lastRow = 0;
		while (cell != null) {
			cell = tracebackRec(X, cell[0], cell[1], alignment, bounds, mode);
			if (cell != null && cell[0] != 0)
				lastCol = cell[1];
			if (cell != null && cell[1] > 3)
				lastRow = cell[0];
		}

		if (borderReached != null)
			return null;

		if (mode == AliMode.FREESHIFT_LEFT) {
			alignment[0] = alignment[0].reverse();
			alignment[1] = alignment[1].reverse();
			alignment[2] = alignment[2].reverse();
		}

		int qRightNotAligned = X[0].length - (startingCell[1] + 3);
		if (mode == AliMode.FREESHIFT_LEFT) {
			int offset = (qRightNotAligned % 3);
			StringBuffer shiftedFrames = new StringBuffer(frames.length());
			for (int i = 0; i < alignment[2].length(); i++) {
				int frame = Integer.parseInt(alignment[2].charAt(i) + "") - 1;
				int shiftedFrame;
				if (frame == offset)
					shiftedFrame = 2 + 1;
				else if (frame == (offset - 1) % 3)
					shiftedFrame = 0 + 1;
				else
					shiftedFrame = 1 + 1;
				shiftedFrames = shiftedFrames.append(shiftedFrame);
			}
			alignment[2] = shiftedFrames;
		}

		int qLeftNotAligned = lastCol < 4 ? 0 : (lastCol - 4);
		int qNotAligned = mode == AliMode.SEMI_GLOBAL ? qLeftNotAligned : qRightNotAligned;
		if (alignment[2].length() > 0)
			qNotAligned -= Integer.parseInt(alignment[2].charAt(0) + "") - 1;
		qNotAligned /= 3;

		int rNotAligned = mode == AliMode.SEMI_GLOBAL ? lastRow - 1 : X.length - (startingCell[0] + 1);
		rNotAligned = rNotAligned < 0 ? 0 : rNotAligned;

		int aliScore = startingCell[2];
		Object[] res = { alignment[1].toString(), alignment[0].toString(), alignment[2].toString(), aliScore, qNotAligned, rNotAligned };

		return res;

	}

	private int[] tracebackRec(int[][] X, int i, int j, StringBuffer[] alignment, Vector<int[]> bounds, AliMode mode) {

		if (verbose) {
			System.out.println(i + " " + j + ": [" + bounds.get(i)[0] + " " + bounds.get(i)[1] + "]");
			System.out.println(alignment[0].toString());
			System.out.println(alignment[1].toString());
			System.out.println(alignment[2].toString());
		}

		// if (bounds != null && j > 4 && (j <= bounds.get(i)[0] || j == bounds.get(i)[1])) {
		if (bounds != null && j > 4 && (j - bounds.get(i)[0] < 2 || j == bounds.get(i)[1])) {
			borderReached = j - bounds.get(i)[0] < 2 ? BORDER_REACHED.LEFT : BORDER_REACHED.RIGHT;
			return null;
		} else if (mode == AliMode.SEMI_GLOBAL && (i == 0 || j < 4))
			return null;

		if (i > 0 || j > 3) {

			int val = X[i][j];

			// checking for match
			if (i > 0 && j > 3) {

				int score = matchScore(i, j, false, mode);

				// checking for hit in same frame
				if ((j - 3) > 0 && X[i - 1][j - 3] == val - score) {
					int[] prevCell = { i - 1, j - 3 };
					alignment[0] = alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1] = alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2] = alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
				}

				// checking for hit in different frames
				if ((j - 2) >= 2 && X[i - 1][j - 2] == val - score + F) {
					int[] prevCell = { i - 1, j - 2 };
					alignment[0] = alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1] = alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2] = alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
				}
				if ((j - 4) >= 4 && X[i - 1][j - 4] == val - score + F) {
					int[] prevCell = { i - 1, j - 4 };
					alignment[0] = alignment[0].insert(0, protein.charAt(i - 1));
					alignment[1] = alignment[1].insert(0, getAA(j - 4, mode));
					alignment[2] = alignment[2].insert(0, ((j - 4) % 3) + 1);
					return prevCell;
				}

			}

			// checking step-wise upper and left direction
			int gap_length = 0;
			while (true) {

				// checking upper direction
				Integer u = null;
				int pos = i - 1 - gap_length;
				if (pos >= 0) {
					if (X[pos][j] - gop - (i - 1 - pos) * gep == val)
						u = pos;
					if (u != null) {
						int[] prevCell = { u, j };
						StringBuffer subSeq = new StringBuffer("");
						for (int k = i; k > u; k--)
							subSeq.insert(0, protein.charAt(k - 1));
						alignment[0] = alignment[0].insert(0, subSeq);
						alignment[1] = alignment[1].insert(0, gapString(i - u));
						if (j > 3)
							alignment[2] = alignment[2].insert(0, frameString(i - u, ((j - 4) % 3) + 1));
						else
							alignment[2] = alignment[2].insert(0, frameString(i - u, j));
						return prevCell;
					}
				}

				// checking left direction
				Integer t = null;
				pos = j - 3 - (gap_length * 3);
				if (pos >= 0) {
					int c = (j - pos) / 3;
					if (X[i][pos] - gop - (c - 1) * gep == val)
						t = pos;
					if (t != null) {
						int[] prevCell = { i, t };
						StringBuffer subSeq = new StringBuffer("");
						for (int k = j; k > t; k -= 3)
							subSeq.insert(0, getAA(k - 4, mode));
						alignment[0] = alignment[0].insert(0, gapString((j - t) / 3));
						alignment[1] = alignment[1].insert(0, subSeq);
						alignment[2] = alignment[2].insert(0, frameString((j - t) / 3, ((j - 4) % 3) + 1));
						return prevCell;
					}
				}

				gap_length++;

			}

		}

		return null;

	}

	private Object gapString(int l) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < l; i++)
			buf.append("-");
		return buf.toString();
	}

	private Object frameString(int l, int frame) {
		StringBuffer buf = new StringBuffer();
		for (int i = 0; i < l; i++)
			buf.append(frame);
		return buf.toString();
	}

	private void initMatrices(int n, int m, int[][] X, int[][] Y, AliMode mode) {

		// initializing everything
		for (int i = 1; i < n; i++) {
			for (int j = 4; j < m - 2; j++) {
				X[i][j] = MIN_VALUE;
			}
		}

		// initializing first column
		for (int i = 1; i < n; i++) {
			X[i][0] = MIN_VALUE;
		}

		// initializing first row
		for (int j = 4; j < m - 2; j++) {
			X[0][j] = mode != AliMode.SEMI_GLOBAL ? -gop - (((j - 1) / 3) - 1) * gep : 0;
		}

		// initializing column 1-3
		X[0][0] = X[0][1] = X[0][2] = X[0][3] = 0;
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < 4; j++) {
				X[i][j] = mode != AliMode.SEMI_GLOBAL ? -gop - ((i - 1) * gep) : 0;
			}
		}

		Y[0][0] = MIN_VALUE;
		Y[0][1] = Y[0][2] = Y[0][3] = 0;
		for (int j = 4; j < m; j++) {
			Y[0][j] = MIN_VALUE;
		}

	}

	private int matchScore(int i, int j, boolean debug, AliMode mode) {
		char a = protein.charAt(i - 1);
		char b = mode != AliMode.FREESHIFT_LEFT ? codonTranslator.translateCodon(dna.charAt(j - 4), dna.charAt(j - 3), dna.charAt(j - 2))
				: codonTranslator.translateCodon(dna.charAt(j - 2), dna.charAt(j - 3), dna.charAt(j - 4));
		return scoringMatrix.getScore(a, b);
	}

	private char getAA(int i, AliMode mode) {

		StringBuffer codon = new StringBuffer(3);
		for (int pos = i; pos < i + 3; pos++)
			codon = codon.append(dna.charAt(pos));
		char aa = mode != AliMode.FREESHIFT_LEFT ? codonTranslator.translateCodon(codon.toString())
				: codonTranslator.translateCodon(codon.reverse().toString());

		return aa;
	}

	// save summation taking care of overflow
	private int add(int a, int b) {
		int res = a + b;
		if (b >= 0 && res < a)
			return MIN_VALUE;
		if (b < 0 && res > a)
			return MIN_VALUE;
		return res;
	}

	// save subtraction taking care of underflow
	private int sub(int a, int b) {
		int res = a - b;
		if (b >= 0 && res > a)
			return MIN_VALUE;
		if (b < 0 && res < a)
			return MIN_VALUE;
		return res;
	}

	// computing most diagonal path
	private Vector<int[]> cmpDiagonal(int n, int m) {

		boolean swapped = false;
		if (m < n) {
			int mem = m;
			m = n;
			n = mem;
			swapped = true;
		}

		Vector<int[]> diagonal = new Vector<int[]>();
		int diff = m - n;
		double frac = diff != 0 ? (double) n / (double) diff : 0.;
		int k = -1;
		if (frac >= 1) {
			int step = (int) Math.round(frac);
			for (int i = 0; i < n; i++) {
				k = k < m - 1 ? k + 1 : k;
				int[] p1 = { i, k };
				diagonal.add(p1);
				if (step != 0 && i % step == 0) {
					k = k < m - 1 ? k + 1 : k;
					int[] p2 = { i, k };
					diagonal.add(p2);
				}
			}
		} else {
			int step = (int) Math.round(1. / frac);
			for (int i = 0; i < n; i++) {
				k = k < m - 1 ? k + 1 : k;
				int[] p1 = { i, k };
				diagonal.add(p1);
				for (int j = 0; j < step; j++) {
					k = k < m - 1 ? k + 1 : k;
					int[] p2 = { i, k };
					diagonal.add(p2);
				}
			}
		}
		for (int j = k + 1; j < m; j++) {
			int[] p = { n - 1, j };
			diagonal.add(p);
		}

		if (swapped) {
			Vector<int[]> realDiagonal = new Vector<int[]>();
			for (int[] p : diagonal) {
				int[] pPrime = { p[1], p[0] };
				realDiagonal.add(pPrime);
			}
			return realDiagonal;
		}

		return diagonal;

	}

	// *********************************************************************

	private void printResult(Object[] res) {

		String ref = (String) res[0];
		String query = (String) res[1];
		String frames = (String) res[2];
		int score = (int) res[3];
		int qAligned = (int) res[4];

		System.out.println(ref + "\n" + query + "\n" + frames + "\n" + score + "\n" + qAligned);

	}

	private void printMatrix(int[][] m, String s1, String s2) {

		System.out.print("\t\t");
		for (int i = 0; i < s2.length() + 4; i++)
			System.out.print(i + "\t");
		System.out.println();

		System.out.print("\t\t\t\t\t\t");
		for (int i = 0; i < s2.length(); i++)
			System.out.print(s2.charAt(i) + "\t");
		System.out.println();

		for (int i = 0; i < m.length; i++) {
			System.out.print(i + "\t");
			if (i > 0)
				System.out.print(s1.charAt(i - 1) + "\t");
			else
				System.out.print("\t");
			for (int j = 0; j < m[0].length; j++) {
				System.out.print(m[i][j] + "\t");
			}
			System.out.println();
		}

	}

}
