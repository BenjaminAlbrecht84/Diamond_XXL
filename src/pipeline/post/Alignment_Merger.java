package pipeline.post;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.BitSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import io.daa.DAA_Reader;
import pipeline.post.Hit.HitType;
import util.CompressAlignment;
import util.DAACompressAlignment;
import util.ReconstructAlignment;
import util.ScoringMatrix;

public class Alignment_Merger {

	private double lambda = 0.267, k = 0.041;
	private File sam_file;
	private DAA_Reader daaReader;

	public Alignment_Merger(HitRun_Rater scorer, File sam_file, DAA_Reader daaReader) {
		this.lambda = scorer.getLambda();
		this.k = scorer.getK();
		this.sam_file = sam_file;
		this.daaReader = daaReader;
	}

	public Hit mergeTwoHits(Hit h1, Hit h2, ScoringMatrix matrix, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {

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
		// int o = r1 - l2 + 1 + h2.numOfQueryInsertions(0, r1 - l2 + 1);
		int o1 = r1 - l2 + 1 + h1.numOfQueryInsertions(l + 1, h1.getAliLength());
		int o2 = r1 - l2 + 1 + h2.numOfQueryInsertions(0, r1 - l2 + 1);
		int r = r2 - r1 + h2.numOfQueryInsertions(o2, h2.getAliLength());
		int length = l + o1 + r;

		// computing new alignmentScores & rawScore
		int rawScore = 0;
		int[] scores = new int[length];
		BitSet query_insertions = new BitSet(length);
		BitSet query_deletions = new BitSet(length);
		for (int i = 0; i < h1.getAliLength(); i++) {
			scores[i] = scores1[i];
			rawScore += scores[i];
			query_insertions.set(i, h1.isQueryInsertion(i));
			query_deletions.set(i, h1.isQueryDeletion(i));
		}

		for (int i = o2; i < h2.getAliLength(); i++) {
			int index = h1.getAliLength() + i - o2;
			scores[index] = scores2[i];
			rawScore += scores[index];
			query_insertions.set(index, h2.isQueryInsertion(i));
			query_deletions.set(index, h2.isQueryDeletion(i));
		}

		int queryLength = h1.getQuery_length() + h2.getQuery_length() - o2;

		// computing new bit score
		double s = rawScore;
		double sPrime = (lambda * s - Math.log(k)) / Math.log(2);
		int bitScore = (int) Math.round(sPrime);

		// generating new alignment strings
		String[] aliStrings1 = daaReader == null ? h1.getAlignmentStrings(rafSAM) : h1.getAlignmentStrings(rafDAA, daaReader);
		String[] aliStrings2 = daaReader == null ? h2.getAlignmentStrings(rafSAM) : h2.getAlignmentStrings(rafDAA, daaReader);
		String[] mergeResult = mergeAliStrings(aliStrings1, aliStrings2, o2, matrix);

		// initializing new hit
		Hit h = new Hit(h1.getId(), h1.getRef_start(), h2.getRef_end(), bitScore, rawScore, h1.getFile_pointer(), h1.getAccessPoint(),
				h1.getQuery_start(), h1.getRef_length(), queryLength, scores, query_insertions, query_deletions, h1.getSubjectID());
		String[] aliStrings = { mergeResult[0], mergeResult[1], mergeResult[2] };
		h.copyAliStrings(aliStrings);
		h.setHitType(HitType.Merged);
		h.setFrame(h1.getFrame());

		Object[] alis = new ReconstructAlignment(matrix).run(mergeResult[1], mergeResult[0], mergeResult[2]);
		String[] ali = { (String) alis[4], (String) alis[5] };
		for (int k = 0; k < ali[0].length(); k++) {
			char c1 = ali[0].charAt(k);
			char c2 = ali[1].charAt(k);
			if (c1 == '-' && c2 == '-') {
				System.out.println(mergeResult[3] + "\n" + mergeResult[4]);
				System.out.println(mergeResult[0]);
				System.out.println(mergeResult[1]);
				System.out.println(mergeResult[2]);
				System.out.println(ali[0] + "\n" + ali[1] + "\n");
			}
		}

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

	private String[] mergeAliStrings(String[] aliStrings1, String[] aliStrings2, int o, ScoringMatrix matrix) {

		// reconstructing both alignments
		Object[] res1 = new ReconstructAlignment(matrix).run(aliStrings1[1], aliStrings1[0], aliStrings1[2]);
		String[] ali1 = { (String) res1[4], (String) res1[5] };
		Object[] res2 = new ReconstructAlignment(matrix).run(aliStrings2[1], aliStrings2[0], aliStrings2[2]);
		String[] ali2 = { (String) res2[4], (String) res2[5] };

		// merging both alignments
		String[] ali = mergeAlignments(ali1, ali2, o);

		// compressing merged alignment
		String[] aliStrings = new CompressAlignment().run(ali);

		String[] result = { aliStrings[0], aliStrings[1], aliStrings[2], ali[0], ali[1] };
		return result;

	}

	private String[] mergeAlignments(String[] ali1, String[] ali2, int o) {
		String s1 = ali1[0].concat(ali2[0].substring(o));
		String s2 = ali1[1].concat(ali2[1].substring(o));
		String[] ali = { s1, s2 };
		return ali;
	}

}
