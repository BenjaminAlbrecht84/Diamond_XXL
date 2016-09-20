package pipeline.post;

import java.io.RandomAccessFile;
import java.util.Vector;

import io.daa.DAA_Reader;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.ScoringMatrix;
import util.SparseString;

public class Hit_Run {

	private boolean isCompleted = false;
	private String readID;
	private Vector<Hit> hitRun;
	private SparseString gi;
	private int sumScore, length, rawScore, coverage, runLength;
	private double sumProbability;
	private Frame_Direction frameDirection;
	private int[] aliStats;

	public Hit_Run(Vector<Hit> hitRun, String readID, SparseString gi, int sumScore, int length, int rawScore, Frame_Direction frameDirection,
			int runLength, Double eValue) {
		this.hitRun = hitRun;
		this.readID = readID;
		this.gi = gi;
		this.sumScore = sumScore;
		this.length = length;
		this.rawScore = rawScore;
		this.frameDirection = frameDirection;
		this.runLength = runLength;
		this.sumProbability = eValue;
		coverage = getCoverage();
	}

	private int getCoverage() {
		double refLength = hitRun.get(0).getRef_length();
		int coverage = (int) Math.floor(((double) runLength / refLength) * 100.);
		return coverage;
	}

	public void update(HitRun_Rater rater, RandomAccessFile rafSAM, RandomAccessFile rafDAA, DAA_Reader daaReader, ScoringMatrix matrix) {
		Object[] res = rater.run(hitRun, frameDirection, rafSAM, rafDAA, readID);
		sumScore = (int) res[0];
		length = (int) res[1];
		rawScore = (int) res[2];
		runLength = (int) res[3];
		sumProbability = (double) res[4];
		coverage = getCoverage();

		aliStats = new int[7];
		for (Hit h : hitRun) {
			int[] stats = h.getAccessPoint() != null ? h.getAlignmentStats(matrix, rafDAA, daaReader) : h.getAlignmentStats(matrix, rafSAM);
			for (int i = 0; i < stats.length; i++)
				aliStats[i] += stats[i];
		}

	}

	public int getCoverge() {
		return coverage;
	}

	public Vector<Hit> getHitRun() {
		return hitRun;
	}

	public SparseString getGi() {
		return gi;
	}

	public int getSumScore() {
		return sumScore;
	}

	public String getReadID() {
		return readID;
	}

	public int getLength() {
		return length;
	}

	public int getRunLength() {
		return runLength;
	}

	public int getRawScore() {
		return rawScore;
	}

	public Frame_Direction getFrameDirection() {
		return frameDirection;
	}

	public boolean isCompleted() {
		return isCompleted;
	}

	public void setCompleted(boolean isCompleted) {
		this.isCompleted = isCompleted;
	}

	public void freeMemory() {
		hitRun.clear();
		hitRun = null;
	}

	public double getEValue() {
		return sumProbability;
	}

	public int[] getAliStats() {
		return aliStats;
	}

}
