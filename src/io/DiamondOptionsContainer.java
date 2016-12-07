package io;

import java.io.File;

public class DiamondOptionsContainer {
	
	final static String[] MATRICIES = { "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM250", "PAM70", "PAM30" };

	// general & IO options
	private int threads = Runtime.getRuntime().availableProcessors();
	private boolean salltitles = false;

	// sensitivity & speed options
	private String[] sensitivity = { "", "sensitive", "more-sensitive" };
	private int sensitivityDegree = 0;

	// scoring & reporting options
	private int gapOpen = 11;
	private int gapExtend = 1;
	private String matrix = "BLOSUM62";
	private boolean seg = true;
	private int maxTargetSeqs = 25;
	private boolean top = false;
	private double eValue = 0.001;
	private Integer minScore = null;
	private Integer queryCover = null;
	private Integer band = null;

	// memory & performance options
	private double blockSize = 2.0;
	private String tmpDir = null;
	private int indexChunks = 4;

	public String getDiamondOptionsString() {
		StringBuffer buf = new StringBuffer();

		// general & IO options
		buf = threads < Runtime.getRuntime().availableProcessors() ? buf.append("-p " + threads + " ") : buf;
		buf = salltitles ? buf.append("--salltitles ") : buf;
		buf = sensitivityDegree > 0 ? buf.append("--" + sensitivity[sensitivityDegree] + " ") : buf;

		// scoring & reporting options
		buf = gapOpen != 11 ? buf.append("--gapopen " + gapOpen + " ") : buf;
		buf = gapExtend != 1 ? buf.append("--gapextend " + gapExtend + " ") : buf;
		buf = !matrix.equals("BLOSUM62") ? buf.append("--matrix " + matrix + " ") : buf;
		buf = !seg ? buf.append("--seq no ") : buf;
		buf = maxTargetSeqs != 25 ? buf.append("--max-target-seqs " + maxTargetSeqs + " ") : buf;
		buf = top ? buf.append("--top ") : buf;
		buf = eValue != 0.001 ? buf.append("--evalue " + eValue + " ") : buf;
		buf = minScore != null ? buf.append("--min-score " + minScore + " ") : buf;
		buf = queryCover != null ? buf.append("--query-cover " + queryCover + " ") : buf;
		buf = band != null ? buf.append("--band " + band + " ") : buf;

		// memory & performance options
		buf = blockSize != 2.0 ? buf.append("--block-size " + blockSize + " ") : buf;
		buf = tmpDir != null ? buf.append("--tmpdir " + tmpDir + " ") : buf;
		buf = indexChunks != 4 ? buf.append("--index-chunks " + indexChunks + " ") : buf;

		return buf.toString();
	}

	public int getThreads() {
		return threads;
	}

	public void setThreads(int threads) {
		this.threads = threads;
	}

	public boolean isSalltitles() {
		return salltitles;
	}

	public void setSalltitles(boolean salltitles) {
		this.salltitles = salltitles;
	}

	public String[] getSensitivity() {
		return sensitivity;
	}

	public void setSensitivity(String[] sensitivity) {
		this.sensitivity = sensitivity;
	}

	public int getSensitivityDegree() {
		return sensitivityDegree;
	}

	public void setSensitivityDegree(int sensitivityDegree) {
		this.sensitivityDegree = sensitivityDegree;
	}

	public int getGapOpen() {
		return gapOpen;
	}

	public void setGapOpen(int gapOpen) {
		this.gapOpen = gapOpen;
	}

	public int getGapExtend() {
		return gapExtend;
	}

	public void setGapExtend(int gapExtend) {
		this.gapExtend = gapExtend;
	}

	public String[] getMatrices() {
		return MATRICIES;
	}

	public String getMatrix() {
		return matrix;
	}

	public boolean setMatrix(String matrix) {
		this.matrix = matrix;
		for (String m : MATRICIES) {
			if (m.equals(matrix))
				return true;
		}
		return false;
	}

	public boolean isSeg() {
		return seg;
	}

	public boolean setSeg(String s) {
		if (s.equals("yes")) {
			seg = true;
			return true;
		}
		if (s.equals("no")) {
			seg = false;
			return true;
		}
		return false;
	}

	public int getMaxTargetSeqs() {
		return maxTargetSeqs;
	}

	public void setMaxTargetSeqs(int maxTargetSeqs) {
		this.maxTargetSeqs = maxTargetSeqs;
	}

	public boolean isTop() {
		return top;
	}

	public void setTop(boolean top) {
		this.top = top;
	}

	public double geteValue() {
		return eValue;
	}

	public void seteValue(double eValue) {
		this.eValue = eValue;
	}

	public Integer getMinScore() {
		return minScore;
	}

	public void setMinScore(Integer minScore) {
		this.minScore = minScore;
	}

	public Integer getQueryCover() {
		return queryCover;
	}

	public void setQueryCover(Integer queryCover) {
		this.queryCover = queryCover;
	}

	public double getBlockSize() {
		return blockSize;
	}

	public void setBlockSize(double blockSize) {
		this.blockSize = blockSize;
	}

	public String getTmpDir() {
		return tmpDir;
	}

	public boolean setTmpDir(String tmpDir) {
		this.tmpDir = tmpDir;
		File tmp = new File(tmpDir);
		if (!tmp.exists() && !tmp.mkdir())
			return false;
		return true;
	}

	public int getIndexChunks() {
		return indexChunks;
	}

	public void setIndexChunks(int indexChunks) {
		this.indexChunks = indexChunks;
	}

	public void setBand(int band) {
		this.band = band;
	}

}
