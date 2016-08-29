package util.frameshiftAligner.banded;

import java.util.Vector;

import util.CodonTranslator;
import util.ScoringMatrix;

public class Banded_Frameshift_Alignment {

	/*
	 * Approach following the work 'Alignments of DNA and protein sequences containing frameshift errors', Guan X. Uberbacher EC., Comput Appl Biosci.
	 * 1996 Feb;12(1):31-40.
	 */

	public enum AliMode {
		GLOBAL, FREESHIFT_LEFT, FREESHIFT_RIGHT, FREESHIFT
	};

	private ScoringMatrix scoringMatrix;
	private int gapOpen, gapExtend, delta;

	private Vector<Integer[]> diagonal;

	public Banded_Frameshift_Alignment(ScoringMatrix scoringMatrix, int delta) {
		this.scoringMatrix = scoringMatrix;
		this.gapOpen = scoringMatrix.getGapOpen();
		this.gapExtend = scoringMatrix.getGapExtend();
		this.delta = delta;
	}

	public Object[] run(String s1, String s2, AliMode mode, double band) {

		String[] frameSeqs = cmpDiffFrameSeqs(s1);
		String s11 = frameSeqs[0];
		String s12 = frameSeqs[1];
		String s13 = frameSeqs[2];

		if (mode == AliMode.FREESHIFT_LEFT) {
			s11 = new StringBuffer(s11).reverse().toString();
			s12 = new StringBuffer(s12).reverse().toString();
			s13 = new StringBuffer(s13).reverse().toString();
			s2 = new StringBuffer(s2).reverse().toString();
		}

		int n1 = s11.length() + 1;
		int n2 = s2.length() + 1;

		// initializing D1, Q1, P1
		int[][] D1 = new int[n1][n2];
		int[][] P1 = new int[n1][n2];
		int[][] Q1 = new int[n1][n2];
		initMatrices(D1, P1, Q1, n1, n2);

		// initializing D2, Q2, P2
		int[][] D2 = new int[n1][n2];
		int[][] P2 = new int[n1][n2];
		int[][] Q2 = new int[n1][n2];
		initMatrices(D2, P2, Q2, n1, n2);

		// initializing D3, Q3, P3
		int[][] D3 = new int[n1][n2];
		int[][] P3 = new int[n1][n2];
		int[][] Q3 = new int[n1][n2];
		initMatrices(D3, P3, Q3, n1, n2);

		// computing diagonal
		if (diagonal == null)
			diagonal = cmpDiagonal(n1, n2);

		// computing bandwidth
		int d = (int) Math.floor((double) Math.max(n1, n2) * band);
		int a = n1 < n2 ? n1 : n2;
		int b = n2 >= n1 ? n2 : n1;
		double frac = (1. - ((double) a / (double) b));
		int dMin = (int) Math.floor((double) Math.max(n1, n2) * frac);
		d += dMin;

		// initializing diagonal pointer
		int pointer = 0;
		while (diagonal.get(pointer)[0] == 0)
			pointer++;

		// initializing diagonal bounds
		Vector<int[]> bounds = new Vector<int[]>();
		int[] initBound = { -1, n2 + 1 };
		bounds.add(initBound);

		// filling matrices
		fillMatrices(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, mode, n1, n2, pointer, d, bounds, s11, s12, s13, s2);

		// performing traceback
		Banded_Frameshift_Traceback traceback = new Banded_Frameshift_Traceback(scoringMatrix, gapOpen, gapExtend, delta);
		int[] startingCell = mode == AliMode.GLOBAL ? getStartingCell_Global(D1, D2, D3, n1, n2) : getStartingCell_FreeShift(D1, D2, D3, n1, n2);
		Object[] alignment = traceback.run(D1, Q1, P1, D2, Q2, P2, D3, Q3, P3, s11, s12, s13, s2, startingCell, bounds, mode);

		// returning result
		if (alignment == null) {
			band = 2. * band <= 1. ? 2. * band : 1.;
			return run(s1, s2, mode, band);
		} else {
			diagonal = null;
			return alignment;
		}

	}

	private void fillMatrices(int[][] D1, int[][] P1, int[][] Q1, int[][] D2, int[][] P2, int[][] Q2, int[][] D3, int[][] P3, int[][] Q3,
			AliMode mode, int n1, int n2, int pointer, int d, Vector<int[]> bounds, String s11, String s12, String s13, String s2) {

		for (int i = 1; i < n1; i++) {

			Integer[] leftCell = diagonal.get(pointer);
			while (pointer < diagonal.size() && diagonal.get(pointer)[0] == i)
				pointer++;
			Integer[] rightCell = diagonal.get(pointer - 1);

			int l = leftCell[1] - d;
			int r = mode == AliMode.GLOBAL ? rightCell[1] + d : n2;
			int[] bound = { l, r };
			bounds.add(bound);

			l = l < 1 ? 1 : l;
			r = r > n2 - 1 ? n2 : r + 1;
			for (int j = l; j < r; j++) {

				int[] res1 = cmpMatrices(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, i, j, s11, s2, 1);
				P1[i][j] = res1[0];
				Q1[i][j] = res1[1];
				D1[i][j] = res1[2];

				int[] res2 = cmpMatrices(D2, P2, Q2, D1, P1, Q1, D3, P3, Q3, i, j, s12, s2, 2);
				P2[i][j] = res2[0];
				Q2[i][j] = res2[1];
				D2[i][j] = res2[2];

				int[] res3 = cmpMatrices(D3, P3, Q3, D1, P1, Q1, D2, P2, Q2, i, j, s13, s2, 3);
				P3[i][j] = res3[0];
				Q3[i][j] = res3[1];
				D3[i][j] = res3[2];

			}

		}

	}

	// dynamic programming approach following Guan&Uberbacher (1996)
	private int[] cmpMatrices(int[][] D1, int[][] P1, int[][] Q1, int[][] D2, int[][] P2, int[][] Q2, int[][] D3, int[][] P3, int[][] Q3, int i,
			int j, String s, String s2, int frame) {

		// boolean debug = (frame == 3 && i == 5 && j == 5);
		// String log = "";

		// checking for staying in frame 1
		int m1 = add(D1[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		P1[i][j] = Math.max(sub(D1[i - 1][j], w(1)), sub(P1[i - 1][j], gapExtend));
		Q1[i][j] = Math.max(sub(D1[i][j - 1], w(1)), sub(Q1[i][j - 1], gapExtend));
		D1[i][j] = Math.max(m1, Math.max(P1[i][j], Q1[i][j]));

		// checking for switching from frame 2 into frame 1
		int m2 = add(D2[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		int t1 = Math.max(sub(D2[i - 1][j], w(1)), sub(P2[i - 1][j], gapExtend));
		int t2 = Math.max(sub(D2[i][j - 1], w(1)), sub(Q2[i][j - 1], gapExtend));
		int t3 = Math.max(m2, Math.max(t1, t2));
		if (sub(t3, delta) > D1[i][j]) {
			P1[i][j] = sub(t1, delta);
			Q1[i][j] = sub(t2, delta);
			D1[i][j] = sub(t3, delta);
		}

		// checking for switching from frame 3 into frame 1
		int m3 = add(D3[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		t1 = Math.max(sub(D3[i - 1][j], w(1)), sub(P3[i - 1][j], gapExtend));
		t2 = Math.max(sub(D3[i][j - 1], w(1)), sub(Q3[i][j - 1], gapExtend));
		t3 = Math.max(m3, Math.max(t1, t2));
		if (sub(t3, delta) > D1[i][j]) {
			P1[i][j] = sub(t1, delta);
			Q1[i][j] = sub(t2, delta);
			D1[i][j] = sub(t3, delta);
		}

		int[] res = { P1[i][j], Q1[i][j], D1[i][j] };
		return res;

	}

	// assessing max cell in bottom right corners
	private int[] getStartingCell_Global(int[][] D1, int[][] D2, int[][] D3, int n1, int n2) {
		int maxF = 1, max = D1[n1 - 1][n2 - 1];
		for (int f = 2; f < 4; f++) {
			int[][] D = f == 2 ? D2 : D3;
			if (D[n1 - 1][n2 - 1] > max) {
				maxF = f;
				max = D[n1 - 1][n2 - 1];
			}
		}
		int[] res = { n1 - 1, n2 - 1, maxF, max };
		return res;
	}

	// assessing max cell in last column
	private int[] getStartingCell_FreeShift(int[][] D1, int[][] D2, int[][] D3, int n1, int n2) {
		int[] max = getMaxInColumn(D1, n2 - 1, 1, n1);
		for (int f = 2; f < 4; f++) {
			int[][] D = f == 2 ? D2 : D3;
			int[] res = getMaxInColumn(D, n2 - 1, f, n1);
			if (res[2] > max[2])
				max = res;
		}
		int[] res = { max[0], max[1], max[3], max[2] };
		return res;
	}

	private int[] getMaxInColumn(int[][] D, int col, int frame, int n1) {
		int[] res = { -1, col, Integer.MIN_VALUE, frame };
		for (int i = 0; i < n1; i++) {
			if (D[i][col] > res[2]) {
				res[0] = i;
				res[2] = D[i][col];
			}
		}
		return res;
	}

	private String[] cmpDiffFrameSeqs(String dna) {

		// compute length of frame sequences
		int l = dna.length(), counter = 0;
		while (counter < 2 || l % 3 != 0) {
			l--;
			counter++;
		}

		// translate DNA sequences
		CodonTranslator translator = new CodonTranslator();
		String s11 = translator.translate(dna.substring(0, l));
		String s12 = translator.translate(dna.substring(1, l + 1));
		String s13 = translator.translate(dna.substring(2, l + 2));

		String[] res = { s11, s12, s13 };
		return res;
	}

	// computing most diagonal path
	private Vector<Integer[]> cmpDiagonal(int n, int m) {

		boolean swapped = false;
		if (m < n) {
			int mem = m;
			m = n;
			n = mem;
			swapped = true;
		}

		Vector<Integer[]> diagonal = new Vector<Integer[]>();
		int diff = m - n;
		double frac = diff != 0 ? (double) n / (double) diff : 0.;
		int k = -1;
		if (frac >= 1) {
			int step = (int) Math.round(frac);
			for (int i = 0; i < n; i++) {
				k = k < m - 1 ? k + 1 : k;
				Integer[] p1 = { i, k };
				diagonal.add(p1);
				if (step != 0 && i % step == 0) {
					k = k < m - 1 ? k + 1 : k;
					Integer[] p2 = { i, k };
					diagonal.add(p2);
				}
			}
		} else {
			int step = (int) Math.round(1. / frac);
			for (int i = 0; i < n; i++) {
				k = k < m - 1 ? k + 1 : k;
				Integer[] p1 = { i, k };
				diagonal.add(p1);
				for (int j = 0; j < step; j++) {
					k = k < m - 1 ? k + 1 : k;
					Integer[] p2 = { i, k };
					diagonal.add(p2);
				}
			}
		}
		for (int j = k + 1; j < m; j++) {
			Integer[] p = { n - 1, j };
			diagonal.add(p);
		}

		if (swapped) {
			Vector<Integer[]> realDiagonal = new Vector<Integer[]>();
			for (Integer[] p : diagonal) {
				Integer[] pPrime = { p[1], p[0] };
				realDiagonal.add(pPrime);
			}
			return realDiagonal;
		}

		return diagonal;

	}

	// standard GOTOH initialization for global alignment
	private void initMatrices(int[][] D, int[][] P, int[][] Q, int n, int m) {

		D[0][0] = 0;
		P[0][0] = Integer.MIN_VALUE;
		Q[0][0] = Integer.MIN_VALUE;

		for (int i = 1; i < n; i++) {
			D[i][0] = -w(i);
			P[i][0] = Integer.MIN_VALUE;
			Q[i][0] = Integer.MIN_VALUE;
		}

		for (int j = 1; j < m; j++) {
			D[0][j] = -w(j);
			P[0][j] = Integer.MIN_VALUE;
			Q[0][j] = Integer.MIN_VALUE;
		}

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < m; j++) {
				D[i][j] = Integer.MIN_VALUE;
				P[i][j] = Integer.MIN_VALUE;
				Q[i][j] = Integer.MIN_VALUE;
			}
		}

	}

	// summing up by taking care of overflow
	private int add(int a, int b) {
		int res = a + b;
		if (b >= 0 && res < a)
			return Integer.MAX_VALUE;
		if (b < 0 && res > a)
			return Integer.MIN_VALUE;
		return res;
	}

	// subtracting by taking care of underflow
	private int sub(int a, int b) {
		int res = a - b;
		if (b >= 0 && res > a)
			return Integer.MIN_VALUE;
		if (b < 0 && res < a)
			return Integer.MAX_VALUE;
		return res;
	}

	private int w(int k) {
		return gapExtend * k + gapOpen;
	}

	// ***********************************************************

	private void printMatrix(Integer[][] M, String s1, String s2) {

		System.out.print("\t\t");
		for (int i = 0; i < s2.length(); i++)
			System.out.print(s2.charAt(i) + "\t");
		System.out.println();

		for (int i = 0; i < M.length; i++) {
			if (i > 0)
				System.out.print(s1.charAt(i - 1) + "\t");
			else
				System.out.print("\t");
			for (int j = 0; j < M[0].length; j++) {
				System.out.print(M[i][j] + "\t");
			}
			System.out.println();
		}

	}

}
