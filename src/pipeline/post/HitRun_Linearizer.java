package pipeline.post;

import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import io.daa.DAA_Reader;
import pipeline.post.mode_one.Alignment_Completer.Alignment_Thread.Hit_Comparator;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.AlignmentEvaluater;
import util.CodonTranslator;
import util.CompressAlignment;
import util.DAACompressAlignment;
import util.ReconstructAlignment;
import util.ScoringMatrix;

public class HitRun_Linearizer {

	private double lambda, k;
	private DAA_Reader daaReader;
	private ScoringMatrix matrix;

	public HitRun_Linearizer(DAA_Reader daaReader, ScoringMatrix matrix, double lambda, double k) {
		this.daaReader = daaReader;
		this.matrix = matrix;
		this.lambda = lambda;
		this.k = k;
	}

	public void run(Hit_Run run, RandomAccessFile rafSAM, RandomAccessFile rafDAA, String q) {

		Vector<Hit> hits = run.getHitRun();

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
		Vector<Hit> linearHits = new Vector<Hit>();
		setHitInfo(h1, rafSAM, rafDAA);
		linearHits.addElement(h1);
		for (int i = 1; i < hits.size(); i++) {

			Hit h2 = hits.get(i);

			int l1 = h1.getRef_start();
			int l2 = h2.getRef_start();
			int r1 = h1.getRef_end();
			int r2 = h2.getRef_end();

			if (r1 <= l2) { // h1 and h2 are consecutive hits

				setHitInfo(h2, rafSAM, rafDAA);
				linearHits.add(h2);
				h1 = h2;

			} else if (r1 >= l2 && r1 < r2) { // h1 and h2 overlap

				int diff = r1 - l2;
				double score1 = h1.getBitScore();
				double score2 = h2.getBitScore();

				if (score1 < score2) {

					// shrinking first hit
					linearHits.removeElement(h1);
					if (h1.getRef_start() < h2.getRef_start()) {

						precomputeAlignmentScore(h1, rafSAM, rafDAA);

						int b = diff + h1.numOfQueryInsertions(h1.getAliLength(), h1.getAliLength() - 1 - diff - 1);
						b = b < h1.getAliLength() ? b : h1.getAliLength() - 1;
						String[] ali1 = getAlignment(h1, rafSAM, rafDAA);

						String[] ali = { ali1[0].substring(0, h1.getAliLength() - 1 - b), ali1[1].substring(0, h1.getAliLength() - 1 - b) };
						int refEnd = h1.getRef_end() - diff - 1;
						int rawScore = matrix.cmpAlignmentScore(ali[0], ali[1]);
						int bitScore = cmpBitScore(rawScore);

						if (ali[0].length() > 0) {

							int[] aliStats = new AlignmentEvaluater().run(ali, matrix);
							int queryLength = aliStats[aliStats.length - 1];
							String[] aliStrings = new CompressAlignment().run(ali);
							Vector<Byte> editOperations = new DAACompressAlignment().run(ali);

							Hit h = new Hit(h1.getId(), h1.getRef_start(), refEnd, bitScore, rawScore, h1.getFile_pointer(), h1.getAccessPoint(),
									h1.getQuery_start(), h1.getRef_length(), queryLength, h1.getSubjectID());
							h.setFrame(h1.getFrame());
							h.setHitType(h1.getHitType());
							Object[] metaInfo = { new Integer(h1.getQuery_start()), new Integer(h1.getFrame()), editOperations };
							h.setMetaInfo(metaInfo);
							h.copyAliStrings(aliStrings);
							h.setAlignmentStats(aliStats);

							// String[] aliString = this.getAlignment(h, rafSAM, rafDAA);
							// System.out.println("1) " + h + "\n" + aliString[0] + "\n" + aliString[1]);

							linearHits.add(h);

						}

					}

					// adopting second hit
					int diff2 = h1.getRef_start() - 1 - h2.getRef_start() > 0 ? h1.getRef_start() - 1 - h2.getRef_start() : 0;
					if (diff2 == 0) {
						setHitInfo(h2, rafSAM, rafDAA);
						linearHits.add(h2);
						h1 = h2;
					} else {

						precomputeAlignmentScore(h2, rafSAM, rafDAA);

						int b = diff2 + 1 + h2.numOfQueryInsertions(0, diff2 + 1);
						String[] ali2 = getAlignment(h2, rafSAM, rafDAA);
						String[] ali = { ali2[0].substring(b), ali2[1].substring(b) };
						int refStart = h2.getRef_start() + diff2 + 1;
						int queryStart = run.getFrameDirection() == Frame_Direction.Positiv
								? h2.getQuery_start() + ((b - h2.numOfQueryDeletionsFixed(0, b)) * 3)
								: h2.getQuery_start() - ((b - h2.numOfQueryDeletionsFixed(0, b)) * 3);

						int rawScore = matrix.cmpAlignmentScore(ali[0], ali[1]);
						int bitScore = cmpBitScore(rawScore);

						int[] aliStats = new AlignmentEvaluater().run(ali, matrix);
						int queryLength = aliStats[aliStats.length - 1];
						String[] aliStrings = new CompressAlignment().run(ali);
						Vector<Byte> editOperations = new DAACompressAlignment().run(ali);

						Hit h = new Hit(h2.getId(), refStart, h2.getRef_end(), bitScore, rawScore, h2.getFile_pointer(), h2.getAccessPoint(),
								queryStart, h2.getRef_length(), queryLength, h2.getSubjectID());
						h.setFrame(h2.getFrame());
						h.setHitType(h2.getHitType());
						Object[] metaInfo = { new Integer(queryStart), new Integer(h2.getFrame()), editOperations };
						h.setMetaInfo(metaInfo);
						h.copyAliStrings(aliStrings);
						h.setAlignmentStats(aliStats);

						linearHits.add(h);

						h1 = h;

						// System.out.println(
						// h.getQuery_start() + " = " + h2.getQuery_start() + " + (" + b + " - " + h2.numOfQueryDeletions(0, diff2) + ")*3");
						// System.out.println(diff2 + " = " + h1.getRef_start() + " -" + 1 + " - " + h2.getRef_start());
						// String[] aliString = this.getAlignment(h, rafSAM, rafDAA);
						// System.out.println("2) " + h + "\n" + aliString[0] + "\n" + aliString[1]);

					}

				} else {

					precomputeAlignmentScore(h2, rafSAM, rafDAA);

					// shrinking second hit
					int b = diff + 1 + h2.numOfQueryInsertions(0, diff + 1);
					String[] ali2 = getAlignment(h2, rafSAM, rafDAA);

					String[] ali = { ali2[0].substring(b), ali2[1].substring(b) };
					int refStart = h2.getRef_start() + diff + 1;
					int queryStart = run.getFrameDirection() == Frame_Direction.Positiv
							? h2.getQuery_start() + ((b - h2.numOfQueryDeletionsFixed(0, b)) * 3)
							: h2.getQuery_start() - ((b - h2.numOfQueryDeletionsFixed(0, b)) * 3);
					int rawScore = matrix.cmpAlignmentScore(ali[0], ali[1]);
					int bitScore = cmpBitScore(rawScore);

					if (ali[0].length() > 0) {

						int[] aliStats = new AlignmentEvaluater().run(ali, matrix);
						int queryLength = aliStats[aliStats.length - 1];
						String[] aliStrings = new CompressAlignment().run(ali);
						Vector<Byte> editOperations = new DAACompressAlignment().run(ali);

						Hit h = new Hit(h2.getId(), refStart, h2.getRef_end(), bitScore, rawScore, h2.getFile_pointer(), h2.getAccessPoint(),
								queryStart, h2.getRef_length(), queryLength, h2.getSubjectID());
						h.setFrame(h2.getFrame());
						h.setHitType(h2.getHitType());
						Object[] metaInfo = { new Integer(queryStart), new Integer(h2.getFrame()), editOperations };
						h.setMetaInfo(metaInfo);
						h.copyAliStrings(aliStrings);
						h.setAlignmentStats(aliStats);

						linearHits.add(h);

						h1 = h;

						// System.out.println(queryStart + " = " + h2.getQuery_start() + " + (" + b + " - " + h2.numOfQueryDeletions(0, b) + ")*3");
						// System.out.println(diff + " = " + r1 + " - " + l2);
						// String[] aliString = this.getAlignment(h, rafSAM, rafDAA);
						// System.out.println("3) " + h + "\n" + aliString[0] + "\n" + aliString[1]);

					}

				}

			}

		}

		// updating hitRun
		Collections.sort(linearHits, new Hit_Comparator());

		run.getHitRun().clear();
		run.getHitRun().addAll(linearHits);

	}

	private void setHitInfo(Hit h, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		if (h.getMetaInfo() == null) {
			String[] ali = getAlignment(h, rafSAM, rafDAA);
			String[] aliStrings = new CompressAlignment().run(ali);
			int[] aliStats = new AlignmentEvaluater().run(ali, matrix);
			Vector<Byte> editOperations = new DAACompressAlignment().run(ali);
			Object[] metaInfo = { h.getQuery_start(), h.getFrame(), editOperations };
			h.setMetaInfo(metaInfo);
			h.copyAliStrings(aliStrings);
			h.setAlignmentStats(aliStats);
		}
	}

	public class Hit_Comparator implements Comparator<Hit> {

		@Override
		public int compare(Hit h1, Hit h2) {
			if (h1.getRef_start() < h2.getRef_start())
				return -1;
			if (h1.getRef_start() > h2.getRef_start())
				return 1;
			return 0;
		}

	}

	private String[] getAlignment(Hit h, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		String[] aliStrings1 = h.getAccessPoint() != null ? h.getAlignmentStrings(rafDAA, daaReader) : h.getAlignmentStrings(rafSAM);
		Object[] res = new ReconstructAlignment(matrix).run(aliStrings1[1], aliStrings1[0], aliStrings1[2]);
		String[] ali = { (String) res[4], (String) res[5] };
		return ali;
	}

	private int cmpBitScore(int rawScore) {
		double s = rawScore;
		double sPrime = (lambda * s - Math.log(k)) / Math.log(2);
		return (int) Math.round(sPrime);
	}

	private void precomputeAlignmentScore(Hit h, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		if (h.getAccessPoint() == null)
			h.getAlignmentScores(matrix, rafSAM);
		else
			h.getAlignmentScores(matrix, rafDAA, daaReader);
	}

}
