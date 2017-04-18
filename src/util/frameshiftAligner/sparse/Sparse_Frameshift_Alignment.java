package util.frameshiftAligner.sparse;

import java.util.Vector;

import util.CodonTranslator;
import util.ScoringMatrix;

public class Sparse_Frameshift_Alignment {

	/*
	 * Approach following the work 'Alignments of DNA and protein sequences containing frameshift errors', Guan X. Uberbacher EC., Comput Appl Biosci.
	 * 1996 Feb;12(1):31-40.
	 */

	private ScoringMatrix scoringMatrix;
	private int gapOpen, gapExtend, delta;

	private int n1, n2;
	private Vector<Integer[]> diagonal;
	private String s1, s2;

	public Sparse_Frameshift_Alignment(ScoringMatrix scoringMatrix, int delta) {
		this.scoringMatrix = scoringMatrix;
		this.gapOpen = scoringMatrix.getGapOpen();
		this.gapExtend = scoringMatrix.getGapExtend();
		this.delta = delta;
	}

	public Object[] run(String s1, String s2, double band) {

		this.s1 = s1;
		this.s2 = s2;

		String[] frameSeqs = cmpDiffFrameSeqs(s1);
		String s11 = frameSeqs[0];
		String s12 = frameSeqs[1];
		String s13 = frameSeqs[2];
		n1 = s11.length() + 1;
		n2 = s2.length() + 1;

		// adjust band
		// int a = n1 <= n2 ? n1 : n2;
		// int b = n2 >= n1 ? n2 : n1;
		// double frac = (1. - ((double) a / (double) b)) / 2.;
		// band = band < frac + band ? frac + band : band;

		// computing diagonal
		if (diagonal == null)
			diagonal = cmpDiagonal(n1, n2);
		int d = (int) Math.ceil((double) Math.min(n1 - 1, n2 - 1) * band);

		// *******************************************

		int n = n1 >= n2 ? n1 : (2 * d) + 1;
		int m = n2 <= n1 ? (2 * d) + 1 : n2;

		// initializing D1, Q1, P1
		Sparse_Matrix D1 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix P1 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix Q1 = new Sparse_Matrix(this, n, m, d);
		initMatrices(D1, P1, Q1, n1, n2);

		// initializing D2, Q2, P2
		Sparse_Matrix D2 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix P2 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix Q2 = new Sparse_Matrix(this, n, m, d);
		initMatrices(D2, P2, Q2, n1, n2);

		// initializing D3, Q3, P3
		Sparse_Matrix D3 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix P3 = new Sparse_Matrix(this, n, m, d);
		Sparse_Matrix Q3 = new Sparse_Matrix(this, n, m, d);
		initMatrices(D3, P3, Q3, n1, n2);

		// *******************************************

		// initializing diagonal bounds
		Vector<int[]> bounds = new Vector<int[]>();
		int[] initBound = { -1, Math.max(n1, n2) + 1 };
		bounds.add(initBound);

		// filling matrices
		if (n1 >= n2)
			fillMatrices_One(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, n1, n2, d, bounds, s11, s12, s13, s2);
		else
			fillMatrices_Two(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, n1, n2, d, bounds, s11, s12, s13, s2);

		// performing traceback
		Sparse_Frameshift_Traceback traceback = new Sparse_Frameshift_Traceback(scoringMatrix, gapOpen, gapExtend, delta, n1, n2);
		int[] startingCell = getStartingCell_Global(D1, D2, D3, n1, n2);
		Object[] alignment = traceback.run(D1, Q1, P1, D2, Q2, P2, D3, Q3, P3, s11, s12, s13, s2, startingCell, bounds);

		// returning result
		if (alignment == null) {
			band = band + 0.1 <= 1. ? band + 0.1 : 1.;
			return run(s1, s2, band);
		} else {
			diagonal = null;
			return alignment;
		}

	}

	private void fillMatrices_One(Sparse_Matrix D1, Sparse_Matrix P1, Sparse_Matrix Q1, Sparse_Matrix D2, Sparse_Matrix P2, Sparse_Matrix Q2,
			Sparse_Matrix D3, Sparse_Matrix P3, Sparse_Matrix Q3, int n1, int n2, int d, Vector<int[]> bounds, String s11, String s12, String s13,
			String s2) {

		int pointer = 1;
		for (int i = 1; i < n1; i++) {

			Integer[] cell = diagonal.get(pointer++);

			int l = cell[1] - d;
			int r = cell[1] + d;
			int[] bound = { l, r };
			bounds.add(bound);

			l = l < 1 ? 1 : l;
			r = r > n2 - 1 ? n2 : r + 1;
			for (int j = l; j < r; j++) {

				int[] res1 = cmpMatrices(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, i, j, s11, s2, 1);
				P1.setValue(i, j, res1[0]);
				Q1.setValue(i, j, res1[1]);
				D1.setValue(i, j, res1[2]);

				int[] res2 = cmpMatrices(D2, P2, Q2, D1, P1, Q1, D3, P3, Q3, i, j, s12, s2, 2);
				P2.setValue(i, j, res2[0]);
				Q2.setValue(i, j, res2[1]);
				D2.setValue(i, j, res2[2]);

				int[] res3 = cmpMatrices(D3, P3, Q3, D1, P1, Q1, D2, P2, Q2, i, j, s13, s2, 3);
				P3.setValue(i, j, res3[0]);
				Q3.setValue(i, j, res3[1]);
				D3.setValue(i, j, res3[2]);

			}

		}

	}

	private void fillMatrices_Two(Sparse_Matrix D1, Sparse_Matrix P1, Sparse_Matrix Q1, Sparse_Matrix D2, Sparse_Matrix P2, Sparse_Matrix Q2,
			Sparse_Matrix D3, Sparse_Matrix P3, Sparse_Matrix Q3, int n1, int n2, int d, Vector<int[]> bounds, String s11, String s12, String s13,
			String s2) {

		int pointer = 1;
		for (int j = 1; j < n2; j++) {

			Integer[] cell = diagonal.get(pointer++);

			int l = cell[0] - d;
			int r = cell[0] + d;
			int[] bound = { l, r };
			bounds.add(bound);

			l = l < 1 ? 1 : l;
			r = r > n1 - 1 ? n1 : r + 1;
			for (int i = l; i < r; i++) {

				int[] res1 = cmpMatrices(D1, P1, Q1, D2, P2, Q2, D3, P3, Q3, i, j, s11, s2, 1);
				P1.setValue(i, j, res1[0]);
				Q1.setValue(i, j, res1[1]);
				D1.setValue(i, j, res1[2]);

				int[] res2 = cmpMatrices(D2, P2, Q2, D1, P1, Q1, D3, P3, Q3, i, j, s12, s2, 2);
				P2.setValue(i, j, res2[0]);
				Q2.setValue(i, j, res2[1]);
				D2.setValue(i, j, res2[2]);

				int[] res3 = cmpMatrices(D3, P3, Q3, D1, P1, Q1, D2, P2, Q2, i, j, s13, s2, 3);
				P3.setValue(i, j, res3[0]);
				Q3.setValue(i, j, res3[1]);
				D3.setValue(i, j, res3[2]);

			}

		}

	}

	// dynamic programming approach following Guan&Uberbacher (1996)
	private int[] cmpMatrices(Sparse_Matrix D1, Sparse_Matrix P1, Sparse_Matrix Q1, Sparse_Matrix D2, Sparse_Matrix P2, Sparse_Matrix Q2,
			Sparse_Matrix D3, Sparse_Matrix P3, Sparse_Matrix Q3, int i, int j, String s, String s2, int frame) {

		// boolean debug = (frame == 3 && i == 5 && j == 5);
		// String log = "";

		// checking for staying in frame 1
		int m1 = add(getEntry(D1.getValue(i - 1, j - 1)), scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		P1.setValue(i, j, Math.max(sub(getEntry(D1.getValue(i - 1, j)), w(1)), sub(getEntry(P1.getValue(i - 1, j)), gapExtend)));
		Q1.setValue(i, j, Math.max(sub(getEntry(D1.getValue(i, j - 1)), w(1)), sub(getEntry(Q1.getValue(i, j - 1)), gapExtend)));
		D1.setValue(i, j, Math.max(m1, Math.max(getEntry(P1.getValue(i, j)), getEntry(Q1.getValue(i, j)))));

		// checking for switching from frame 2 into frame 1
		int m2 = add(getEntry(D2.getValue(i - 1, j - 1)), scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		int t1 = Math.max(sub(getEntry(D2.getValue(i - 1, j)), w(1)), sub(getEntry(P2.getValue(i - 1, j)), gapExtend));
		int t2 = Math.max(sub(getEntry(D2.getValue(i, j - 1)), w(1)), sub(getEntry(Q2.getValue(i, j - 1)), gapExtend));
		int t3 = Math.max(m2, Math.max(t1, t2));
		if (sub(t3, delta) > getEntry(D1.getValue(i, j))) {
			P1.setValue(i, j, sub(t1, delta));
			Q1.setValue(i, j, sub(t2, delta));
			D1.setValue(i, j, sub(t3, delta));
		}

		// checking for switching from frame 3 into frame 1
		int m3 = add(getEntry(D3.getValue(i - 1, j - 1)), scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		t1 = Math.max(sub(getEntry(D3.getValue(i - 1, j)), w(1)), sub(getEntry(P3.getValue(i - 1, j)), gapExtend));
		t2 = Math.max(sub(getEntry(D3.getValue(i, j - 1)), w(1)), sub(getEntry(Q3.getValue(i, j - 1)), gapExtend));
		t3 = Math.max(m3, Math.max(t1, t2));
		if (sub(t3, delta) > getEntry(D1.getValue(i, j))) {
			P1.setValue(i, j, sub(t1, delta));
			Q1.setValue(i, j, sub(t2, delta));
			D1.setValue(i, j, sub(t3, delta));
		}

		int[] res = { getEntry(P1.getValue(i, j)), getEntry(Q1.getValue(i, j)), getEntry(D1.getValue(i, j)) };
		return res;

	}

	private int getEntry(Integer i) {
		if (i != null)
			return i;
		return Integer.MIN_VALUE;
	}

	// assessing max cell in bottom right corners
	private int[] getStartingCell_Global(Sparse_Matrix D1, Sparse_Matrix D2, Sparse_Matrix D3, int n1, int n2) {

		int maxF = 1, max = D1.getValue(n1 - 1, n2 - 1);
		for (int f = 2; f < 4; f++) {
			Sparse_Matrix D = f == 2 ? D2 : D3;
			if (D.getValue(n1 - 1, n2 - 1) > max) {
				maxF = f;
				max = D.getValue(n1 - 1, n2 - 1);
			}
		}
		int[] res = { n1 - 1, n2 - 1, maxF, max };
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
		if (frac >= 1 || frac == 0) {
			int step = (int) Math.ceil(frac);
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
			int step = (int) Math.ceil(1. / frac);
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
	private void initMatrices(Sparse_Matrix D, Sparse_Matrix P, Sparse_Matrix Q, int n, int m) {

		D.setValue(0, 0, 0);
		P.setValue(0, 0, Integer.MIN_VALUE);
		Q.setValue(0, 0, Integer.MIN_VALUE);

		for (int i = 1; i < n; i++) {
			if (!D.checkIndices(i, 0))
				break;
			D.setValue(i, 0, -w(i));
			P.setValue(i, 0, Integer.MIN_VALUE);
			Q.setValue(i, 0, Integer.MIN_VALUE);
		}

		for (int j = 1; j < m; j++) {
			if (!D.checkIndices(0, j))
				break;
			D.setValue(0, j, -w(j));
			P.setValue(0, j, Integer.MIN_VALUE);
			Q.setValue(0, j, Integer.MIN_VALUE);
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

	public int getN1() {
		return n1;
	}

	public int getN2() {
		return n2;
	}

	public Vector<Integer[]> getDiagonal() {
		return diagonal;
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
