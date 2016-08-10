package util.frameshiftAligner;

import util.CodonTranslator;
import util.ScoringMatrix;

public class Frameshift_Alignment {

	/*
	 * Approach following the work 'Alignments of DNA and protein sequences
	 * containing frameshift errors', Guan X. Uberbacher EC., Comput Appl
	 * Biosci. 1996 Feb;12(1):31-40.
	 */

	public enum AliMode {
		GLOBAL, FREESHIFT
	};

	private ScoringMatrix scoringMatrix;
	private int gapOpen, gapExtend, delta;

	public Frameshift_Alignment(ScoringMatrix scoringMatrix, int delta) {
		this.scoringMatrix = scoringMatrix;
		this.gapOpen = scoringMatrix.getGapOpen();
		this.gapExtend = scoringMatrix.getGapExtend();
		this.delta = delta;
	}

	public String[] run(String s1, String s2, AliMode mode) {
		
		String[] frameSeqs = cmpDiffFrameSeqs(s1);
		String s11 = frameSeqs[0];
		String s12 = frameSeqs[1];
		String s13 = frameSeqs[2];

		// initialization of D1, Q1, P1
		int n1 = s11.length() + 1;
		int n2 = s2.length() + 1;
		int[][] D1 = new int[n1][n2];
		int[][] P1 = new int[n1][n2];
		int[][] Q1 = new int[n1][n2];
		initMatrices(D1, P1, Q1, n1, n2);

		// initialization of D2, Q2, P2
		int[][] D2 = new int[n1][n2];
		int[][] P2 = new int[n1][n2];
		int[][] Q2 = new int[n1][n2];
		initMatrices(D2, P2, Q2, n1, n2);

		// initialization of D3, Q3, P3
		int[][] D3 = new int[n1][n2];
		int[][] P3 = new int[n1][n2];
		int[][] Q3 = new int[n1][n2];
		initMatrices(D3, P3, Q3, n1, n2);

		// filling matrices
		for (int i = 1; i < n1; i++) {
			for (int j = 1; j < n2; j++) {

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

		// performing traceback
		Frameshift_Traceback traceback = new Frameshift_Traceback(scoringMatrix, gapOpen, gapExtend, delta);
		int[] startingCell = mode == AliMode.GLOBAL ? getStartingCell_Global(D1, D2, D3, n1, n2)
				: getStartingCell_FreeShift(D1, D2, D3, n2);
		String[] alignment = traceback.run(D1, Q1, P1, D2, Q2, P2, D3, Q3, P3, s11, s12, s13, s2, startingCell);

//		System.out.println("\n" + alignment[0] + "\n" + alignment[1] + "\n" + alignment[2] + " \n" + alignment[3]);

		return alignment;

	}

	// determine max cell in bottom right corners
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

	// determine max cell in last rows
	private int[] getStartingCell_FreeShift(int[][] D1, int[][] D2, int[][] D3, int n2) {
		int[] max = getMaxInColumn(D1, n2 - 1, 1);
		for (int f = 2; f < 4; f++) {
			int[][] D = f == 2 ? D2 : D3;
			int[] res = getMaxInColumn(D, n2 - 1, f);
			if (res[2] > max[2])
				max = res;
		}
		int[] res = { max[0], max[1], max[3], max[2] };
		return res;
	}

	private int[] getMaxInColumn(int[][] D, int col, int frame) {
		int[] res = { -1, col, Integer.MIN_VALUE, frame };
		for (int i = 0; i < D.length; i++) {
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

	// dynamic programming approach following Guan&Uberbacher (1996)
	private int[] cmpMatrices(int[][] D1, int[][] P1, int[][] Q1, int[][] D2, int[][] P2, int[][] Q2, int[][] D3,
			int[][] P3, int[][] Q3, int i, int j, String s, String s2, int frame) {

		// boolean debug = (frame == 3 && i == 5 && j == 5);
		// String log = "";

		// checking for staying in frame 1
		int m1 = add(D1[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		P1[i][j] = Math.max(sub(D1[i - 1][j], w(1)), sub(P1[i - 1][j], gapExtend));
		Q1[i][j] = Math.max(sub(D1[i][j - 1], w(1)), sub(Q1[i][j - 1], gapExtend));
		D1[i][j] = Math.max(m1, Math.max(P1[i][j], Q1[i][j]));

		// checking for shifting from frame 2 into frame 1
		int m2 = add(D2[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		int t1 = Math.max(sub(D2[i - 1][j], w(1)), sub(P2[i - 1][j], gapExtend));
		int t2 = Math.max(sub(D2[i][j - 1], w(1)), sub(Q2[i][j - 1], gapExtend));
		int t3 = Math.max(m2, Math.max(t1, t2));
		if (t3 - delta > D1[i][j]) {
			P1[i][j] = t1 - delta;
			Q1[i][j] = t2 - delta;
			D1[i][j] = t3 - delta;
		}

		// checking for shifting from frame 3 into frame 1
		int m3 = add(D3[i - 1][j - 1], scoringMatrix.getScore(s.charAt(i - 1), s2.charAt(j - 1)));
		t1 = Math.max(sub(D3[i - 1][j], w(1)), sub(P3[i - 1][j], gapExtend));
		t2 = Math.max(sub(D3[i][j - 1], w(1)), sub(Q3[i][j - 1], gapExtend));
		t3 = Math.max(m3, Math.max(t1, t2));
		if (t3 - delta > D1[i][j]) {
			P1[i][j] = t1 - delta;
			Q1[i][j] = t2 - delta;
			D1[i][j] = t3 - delta;
		}

		int[] res = { P1[i][j], Q1[i][j], D1[i][j] };
		return res;

	}

	// basic gotoh initialization for global alignment
	private void initMatrices(int[][] D, int[][] P, int[][] Q, int n1, int n2) {
		D[0][0] = 0;
		P[0][0] = Integer.MIN_VALUE;
		Q[0][0] = Integer.MIN_VALUE;
		for (int i = 1; i < n1; i++) {
			D[i][0] = -w(i);
			P[i][0] = Integer.MIN_VALUE;
			Q[i][0] = Integer.MIN_VALUE;
		}
		for (int j = 1; j < n2; j++) {
			D[0][j] = -w(j);
			P[0][j] = Integer.MIN_VALUE;
			Q[0][j] = Integer.MIN_VALUE;
		}
	}

	// save summation taking care of overflow
	private int add(int a, int b) {
		int res = a + b;
		if (b >= 0 && res < a)
			return Integer.MAX_VALUE;
		if (b < 0 && res > a)
			return Integer.MIN_VALUE;
		return res;
	}

	// save subtraction taking care of underflow
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

	private void printMatrix(int[][] m, String s1, String s2) {

		System.out.print("\t\t");
		for (int i = 0; i < s2.length(); i++)
			System.out.print(s2.charAt(i) + "\t");
		System.out.println();

		for (int i = 0; i < m.length; i++) {
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
