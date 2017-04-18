package pipeline.post;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import io.daa.DAA_Reader;
import pipeline.post.Hit.HitType;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.CompressAlignment;
import util.DAACompressAlignment;
import util.ReconstructAlignment;
import util.ScoringMatrix;
import util.SparseString;

public class Alignment_Merger {

	private double lambda = 0.267, k = 0.041;
	private int step, length;
	private File sam_file;
	private DAA_Reader daaReader;

	public Alignment_Merger(HitRun_Rater scorer, File sam_file, DAA_Reader daaReader, int step, int length) {
		this.lambda = scorer.getLambda();
		this.k = scorer.getK();
		this.sam_file = sam_file;
		this.daaReader = daaReader;
		this.step = step;
		this.length = length;
	}

	public Hit mergeTwoHits(Hit h1, Hit h2, ScoringMatrix matrix, RandomAccessFile rafSAM, RandomAccessFile rafDAA, Hit_Run run,
			HashMap<String, String> readIDToSeq, HashMap<SparseString, String> giToSeq, Frame_Direction dir) {

		if (h1.getRef_start() > h2.getRef_start()) {
			Hit hMem = h1;
			h1 = h2;
			h2 = hMem;
		}

		if (h1.getRef_end() < h2.getRef_start() - 1)
			throw new IllegalArgumentException("Hits have no overlap nor are consecutive!");
		else if (h1.getRef_end() >= h2.getRef_end()) {
			return h1;
		}

		int l1 = h1.getRef_start();
		int r1 = h1.getRef_end();
		int l2 = h2.getRef_start();
		int r2 = h2.getRef_end();

		int[] scores1 = daaReader == null ? h1.getAlignmentScores(matrix, rafSAM) : h1.getAlignmentScores(matrix, rafDAA, daaReader);
		int[] scores2 = daaReader == null ? h2.getAlignmentScores(matrix, rafSAM) : h2.getAlignmentScores(matrix, rafDAA, daaReader);

		int l = l2 - l1 + h1.numOfQueryInsertions(0, l2 - l1 + 1);
		int o1 = r1 - l2 + 1 + h1.numOfQueryInsertions(l + 1, h1.getAliLength());
		int o2 = r1 - l2 + 1 + h2.numOfQueryInsertions(0, r1 - l2 + 1);
		int r = r2 - r1 + h2.numOfQueryInsertions(o2, h2.getAliLength());
		// int length = l + o1 + r;

		// compute cutPoints
		String[] aliStrings1 = daaReader == null ? h1.getAlignmentStrings(rafSAM) : h1.getAlignmentStrings(rafDAA, daaReader);
		String[] aliStrings2 = daaReader == null ? h2.getAlignmentStrings(rafSAM) : h2.getAlignmentStrings(rafDAA, daaReader);
		int[] cutPoints = cmpCutPoint(o2, h1, h2, aliStrings1, aliStrings2, dir, matrix);

		if (cutPoints != null) {
			int length = cutPoints[0] + h2.getAliLength() - cutPoints[1];

			// computing new alignment properties (scores, insertion, deletions,..)
			int rawScore = 0, queryLength = 0;
			int[] scores = new int[length];
			BitSet query_insertions = new BitSet(length);
			BitSet query_deletions = new BitSet(length);
			// for (int i = 0; i < h1.getAliLength(); i++) {
			for (int i = 0; i < cutPoints[0]; i++) {
				scores[i] = scores1[i];
				rawScore += scores[i];
				query_insertions.set(i, h1.isQueryInsertion(i));
				query_deletions.set(i, h1.isQueryDeletion(i));
				if (!h1.isQueryDeletion(i))
					queryLength++;
			}
			// for (int i = o2; i < h2.getAliLength(); i++) {
			for (int i = cutPoints[1]; i < h2.getAliLength(); i++) {
				// int index = h1.getAliLength() + i - o2;
				int index = cutPoints[0] + i - cutPoints[1];
				scores[index] = scores2[i];
				rawScore += scores[index];
				query_insertions.set(index, h2.isQueryInsertion(i));
				query_deletions.set(index, h2.isQueryDeletion(i));
				if (!h2.isQueryDeletion(i))
					queryLength++;
			}
			// int queryLength = h1.getQuery_length() + h2.getQuery_length() - o2;

			// computing new bit score
			double s = rawScore;
			double sPrime = (lambda * s - Math.log(k)) / Math.log(2);
			int bitScore = (int) Math.round(sPrime);

			// generating new alignment strings
			String[] mergeResult = mergeAliStrings(aliStrings1, aliStrings2, o2, matrix, cutPoints);

			// initializing new hit
			Hit h = new Hit(h1.getId(), h1.getRef_start(), h2.getRef_end(), bitScore, rawScore, h1.getFile_pointer(), h1.getAccessPoint(),
					h1.getQuery_start(), h1.getRef_length(), queryLength, scores, query_insertions, query_deletions, h1.getSubjectID());
			String[] aliStrings = { mergeResult[0], mergeResult[1], mergeResult[2] };
			h.copyAliStrings(aliStrings);
			h.setHitType(HitType.Merged);
			h.setFrame(h1.getFrame());

			if (h1.getMetaInfo() != null) {
				Object[] meta = new Object[3];
				meta[0] = new Integer((int) h1.getMetaInfo()[0]);
				meta[1] = new Integer((int) h1.getMetaInfo()[1]);
				String[] mergedAli = { mergeResult[3], mergeResult[4] };
				meta[2] = new DAACompressAlignment().run(mergedAli);
				h.setMetaInfo(meta);
			}

			return h;

		}

		return null;

	}

	private int[] cmpCutPoint(int o2, Hit h1, Hit h2, String[] aliStrings1, String[] aliStrings2, Frame_Direction dir, ScoringMatrix matrix) {

		// reconstructing both alignments
		Object[] res1 = new ReconstructAlignment(matrix).run(aliStrings1[1], aliStrings1[0], aliStrings1[2]);
		String[] ali1 = { (String) res1[4], (String) res1[5] };
		Object[] res2 = new ReconstructAlignment(matrix).run(aliStrings2[1], aliStrings2[0], aliStrings2[2]);
		String[] ali2 = { (String) res2[4], (String) res2[5] };

		// initializing putative cutting points
		int c1 = ali1[0].length();
		int c2 = o2;

		// refining cutting points
		int qStart1 = cmpQueryStart(h1, dir);
		int qStart2 = cmpQueryStart(h2, dir);

		if (c1 > 0 && c2 > 0) {
			int q1 = ali1[0].charAt(c1 - 1) == '-' ? -1 : qStart1 + ((c1) - h1.numOfQueryDeletionsFixed(0, c1)) * 3;
			int q2 = ali2[0].charAt(c2 - 1) == '-' ? -1 : qStart2 + ((c2) - h2.numOfQueryDeletionsFixed(0, c2)) * 3;
			while (c1 > 1 && c2 > 1 && (q1 != q2 || q1 == -1 || q2 == -1)) {
				if (ali2[1].charAt(c2 - 1) != '-')
					c1--;
				if (ali1[1].charAt(c1 - 1) != '-')
					c2--;
				q1 = ali1[0].charAt(c1 - 1) == '-' ? -1 : qStart1 + ((c1) - h1.numOfQueryDeletionsFixed(0, c1)) * 3;
				q2 = ali2[0].charAt(c2 - 1) == '-' ? -1 : qStart2 + ((c2) - h2.numOfQueryDeletionsFixed(0, c2)) * 3;
			}
			if (q1 != q2)
				return null;
		}

		int[] cutPoints = { c1, c2 };
		return cutPoints;
	}

	private int cmpQueryStart(Hit h, Frame_Direction dir) {
		int offset = Math.abs(h.getId() * step);
		int qStart = dir == Frame_Direction.Positiv ? h.getQuery_start() + offset - 1 : length - h.getQuery_start() - offset;
		return qStart;
	}

	private String[] mergeAliStrings(String[] aliStrings1, String[] aliStrings2, int o, ScoringMatrix matrix, int[] cutPoints) {

		// reconstructing both alignments
		Object[] res1 = new ReconstructAlignment(matrix).run(aliStrings1[1], aliStrings1[0], aliStrings1[2]);
		String[] ali1 = { (String) res1[4], (String) res1[5] };
		Object[] res2 = new ReconstructAlignment(matrix).run(aliStrings2[1], aliStrings2[0], aliStrings2[2]);
		String[] ali2 = { (String) res2[4], (String) res2[5] };

		// merging both alignments
		// String[] ali = mergeAlignments(ali1, ali2, o);
		String[] ali = mergeAlignments(ali1, ali2, cutPoints);

		// compressing merged alignment
		String[] aliStrings = new CompressAlignment().run(ali);

		String[] result = { aliStrings[0], aliStrings[1], aliStrings[2], ali[0], ali[1] };
		return result;

	}

	private String[] mergeAlignments(String[] ali1, String[] ali2, int[] cutPoints) {
		String s1 = ali1[0].substring(0, cutPoints[0]).concat(ali2[0].substring(cutPoints[1]));
		String s2 = ali1[1].substring(0, cutPoints[0]).concat(ali2[1].substring(cutPoints[1]));
		String[] ali = { s1, s2 };
		return ali;
	}

	private String[] mergeAlignments(String[] ali1, String[] ali2, int o) {
		String s1 = ali1[0].concat(ali2[0].substring(o));
		String s2 = ali1[1].concat(ali2[1].substring(o));
		String[] ali = { s1, s2 };
		return ali;
	}

}
