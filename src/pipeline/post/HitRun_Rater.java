package pipeline.post;

import java.io.RandomAccessFile;
import java.math.BigInteger;
import java.util.HashMap;
import java.util.Vector;

import io.daa.DAA_Reader;
import pipeline.post.Hit.HitType;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.ScoringMatrix;

public class HitRun_Rater {

	// gap-decay constant of BlastX
	private final double beta = 0.1;

	private double lambda, k;
	private BigInteger m;
	private DAA_Reader daaReader;
	private ScoringMatrix matrix;
	private HashMap<String, Integer> readToLength;

	public HitRun_Rater(double lambda, double k, BigInteger m, DAA_Reader daaReader, ScoringMatrix matrix, HashMap<String, Integer> readToLength) {
		this.lambda = lambda;
		this.k = k;
		this.m = m;
		this.daaReader = daaReader;
		this.matrix = matrix;
		this.readToLength = readToLength;
	}

	public Object[] run(Vector<Hit> hits, Frame_Direction dir, RandomAccessFile rafSAM, RandomAccessFile rafDAA, String readID, int gi) {

		if (hits.isEmpty()) {
			Object[] result = { 0, 0, 0, 0, 0 };
			return result;
		}

		if (hits.size() == 1) {
			Hit h = hits.firstElement();
			double sPrime = h.getBitScore();
			double n = readToLength != null ? readToLength.get(readID) : hits.get(0).getRef_length();
			double eValue = m.doubleValue() * n * Math.pow(2, -sPrime);
			Object[] result = { h.getBitScore(), h.getQuery_length(), h.getRawScore(), refCoverage(h), eValue };
			return result;
		}

		// run of hits has to be in increasing order
		for (int i = 1; i < hits.size(); i++) {
			Hit h1 = hits.get(i - 1);
			Hit h2 = hits.get(i);
			if (h1.getRef_start() > h2.getRef_start()) {
				throw new IllegalArgumentException("ERROR: Run of hits is not in increasing order!");
			}
		}

		// determining properties of run
		Hit h1 = hits.get(0);
		double hit_length = h1.getQuery_length(), sum = h1.getRawScore(), g = 1, r = hits.size();
		int refCover = refCoverage(h1);
		for (int i = 1; i < hits.size(); i++) {

			Hit h2 = hits.get(i);

			int l1 = h1.getRef_start();
			int l2 = h2.getRef_start();
			int r1 = h1.getRef_end();
			int r2 = h2.getRef_end();

			if (r1 <= l2) { // h1 and h2 are consecutive hits
				refCover += refCoverage(h2);
				hit_length += h2.getQuery_length();
				sum += h2.getRawScore();
				g += l2 - r1 + 1;
				h1 = h2;
			} else if (r1 >= l2 && r1 < r2) { // h1 and h2 overlap

				precomputeAlignmentScore(h1, rafSAM, rafDAA);
				precomputeAlignmentScore(h2, rafSAM, rafDAA);

				int diff = r1 - l2;
				// int diff = r1 - l2 + 1;
				refCover += refCoverage(h2) - diff;

				double score1 = h1.getBitScore();
				double score2 = h2.getBitScore();

				if (score1 < score2) {
					int b = diff + h1.numOfQueryInsertions(r1 - diff, r1);
					b = b < h1.getAlignmentScores(matrix, null).length ? b : h1.getAlignmentScores(matrix, null).length - 1;
					hit_length += h2.getQuery_length() - b;
					for (int j = 0; j < b; j++) 
						sum -= h1.getAlignmentScores(matrix, null)[h1.getAlignmentScores(matrix, null).length - 1 - j];
					for (int j = 0; j < h2.getAlignmentScores(matrix, null).length; j++)
						sum += h2.getAlignmentScores(matrix, null)[j];
				} else {
					int b = diff + h2.numOfQueryInsertions(0, diff);
					hit_length += h2.getQuery_length() - b;
					for (int j = b; j < h2.getAlignmentScores(matrix, null).length; j++)
						sum += h2.getAlignmentScores(matrix, null)[j];
				}

				h1 = h2;

			} else { // h2 is contained in h1
				r--;
			}

		}

		double refLength = hits.get(0).getRef_length();
		int coverage = (int) Math.floor(((double) refCover / refLength) * 100.);

		if (coverage < 95) {

			// applying Equation 4-17/4-18 given in the Blast-Book of o'Reilly
			int sumScore = -1;
			double n = readToLength != null ? readToLength.get(readID) : hits.get(0).getRef_length();
			if (g != 0)
				sumScore = (int) Math.round(
						lambda * sum - Math.log(k * m.doubleValue() * n) - (r - 1) * (Math.log(k) + 2 * Math.log(g)) - Math.log10(factorial(r)));
			else
				sumScore = (int) Math.round(lambda * sum - r * Math.log(k * m.doubleValue() * n) + Math.log(factorial(r)));

			// applying Equation 4.19-4.21 given in the Blast-Book of o'Reilly
			double P_r = (Math.pow(Math.E, -(double) sumScore) * Math.pow((double) sumScore, r - 1)) / (factorial(r) * factorial(r - 1));
			double P_r_prime = P_r / Math.pow(beta, r - 1) * (1 - beta);
			double eValue = (m.doubleValue() / n) * P_r_prime;

			Object[] result = { sumScore, (int) hit_length, (int) sum, refCover, eValue };
			return result;

		} else {

			double rawScore = sum;

			double sPrime = (lambda * rawScore - Math.log(k)) / Math.log(2);
			int bitScore = (int) Math.round(sPrime);

			double n = readToLength != null ? readToLength.get(readID) : hits.get(0).getRef_length();
			double eValue = m.doubleValue() * n * Math.pow(2, -sPrime);

			Object[] result = { bitScore, (int) hit_length, (int) rawScore, refCover, eValue };
			return result;

		}

	}

	private void precomputeAlignmentScore(Hit h, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		if (h.getAccessPoint() == null)
			h.getAlignmentScores(matrix, rafSAM);
		else
			h.getAlignmentScores(matrix, rafDAA, daaReader);
	}

	private int refCoverage(Hit h) {
		return h.getRef_end() - h.getRef_start();
	}

	private double factorial(double n) {
		double fact = 1.;
		for (int i = 1; i <= n; i++)
			fact *= (double) i;
		return (double) fact;
	}

	public double getLambda() {
		return lambda;
	}

	public double getK() {
		return k;
	}

	public BigInteger getM() {
		return m;
	}

}
