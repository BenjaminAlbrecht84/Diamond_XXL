package util.frameshiftAligner.sparse;

public class Sparse_Matrix {

	private Sparse_Frameshift_Alignment ali;

	private Integer[][] M;
	private int mid;

	public Sparse_Matrix(Sparse_Frameshift_Alignment ali, int n, int m, int mid) {
		this.ali = ali;
		M = new Integer[n][m];
		// mid = n >= m ? (int) Math.floor((double) m / 2.) : (int) Math.floor((double) n / 2.);
		this.mid = mid;
	}

	public Integer getValue(int i, int j) {
		int[] m = mapIndices(i, j);
		if (rangeCheck(m))
			return M[m[0]][m[1]];
		return null;
	}

	public void setValue(int i, int j, int value) {
		int[] m = mapIndices(i, j);
		if (rangeCheck(m))
			M[m[0]][m[1]] = value;
	}

	public boolean checkIndices(int i, int j) {
		int[] m = mapIndices(i, j);
		if (m[0] >= 0 && m[0] < M.length && m[1] >= 0 && m[1] < M[0].length)
			return true;
		return false;
	}

	private boolean rangeCheck(int[] m) {
		if (m[0] >= 0 && m[0] < M.length && m[1] >= 0 && m[1] < M[0].length)
			return true;
		return false;
	}

	public int[] mapIndices(int i, int j) {
		if (ali.getN1() < ali.getN2()) {
			int mid_i = ali.getDiagonal().get(j)[0];
			int i_prime = mid + i - mid_i;
			int[] mappedIndices = { i_prime, j };
			return mappedIndices;
		} else {
			int mid_j = ali.getDiagonal().get(i)[1];
			int j_prime = mid + j - mid_j;
			int[] mappedIndices = { i, j_prime };
			return mappedIndices;
		}
	}

	public Integer[][] getMatrix() {
		return M;
	}

	public void printMatrix() {
		System.out.println();
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[0].length; j++) {
				System.out.print(M[i][j] + "\t");
			}
			System.out.println();
		}
	}

}
