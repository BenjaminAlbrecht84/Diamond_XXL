package pipeline.post;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.BitSet;

import io.daa.DAA_Reader;
import pipeline.post.Hit.HitType;
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

		// System.out.println(h1.getRef_start() + " " + h2.getRef_start());

		if (h1.getRef_end() < h2.getRef_start())
			throw new IllegalArgumentException("Hits have no overlap!");
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

		int query_length = h1.getQuery_length() + h2.getQuery_length() - o2 + h2.numOfQueryDeletions(0, o2);

		// computing new bit score
		double s = rawScore;
		double sPrime = (lambda * s - Math.log(k)) / Math.log(2);
		int bitScore = (int) Math.round(sPrime);

		// generating new alignment strings
		String[] aliStrings1 = daaReader == null ? h1.getAlignmentStrings(rafSAM) : h1.getAlignmentStrings(rafDAA, daaReader);
		String[] aliStrings2 = daaReader == null ? h2.getAlignmentStrings(rafSAM) : h2.getAlignmentStrings(rafDAA, daaReader);
		String[] aliStrings = mergeAliStrings(aliStrings1, aliStrings2, o2, matrix);

		// initialize new hit
		Hit h = new Hit(h1.getId(), h1.getRef_start(), h2.getRef_end(), bitScore, rawScore, h1.getFile_pointer(), h1.getAccessPoint(),
				h1.getQuery_start(), h1.getRef_length(), query_length, scores, query_insertions, query_deletions);
		h.copyAliStrings(aliStrings);
		h.setHitType(HitType.Merged);

		return h;

	}

	private String[] mergeAliStrings(String[] aliStrings1, String[] aliStrings2, int o, ScoringMatrix matrix) {

		// reconstruct both alignments
		Object[] res1 = new ReconstructAlignment(matrix).run(aliStrings1[1], aliStrings1[0], aliStrings1[2]);
		String[] ali1 = { (String) res1[4], (String) res1[5] };
		Object[] res2 = new ReconstructAlignment(matrix).run(aliStrings2[1], aliStrings2[0], aliStrings2[2]);
		String[] ali2 = { (String) res2[4], (String) res2[5] };

		// merge both alignments
		String[] ali = mergeAlignments(ali1, ali2, o);

		// compute alignment strings
		String[] aliStrings = new String[3];
		aliStrings[0] = deriveCIGAR(ali);
		aliStrings[1] = deriveSEQ(ali);
		aliStrings[2] = deriveMD(ali);

		// System.out.println(ali[0] + "\n" + ali[1]);
		// System.out.println(aliStrings[0] + "\n" + aliStrings[1] + "\n" +
		// aliStrings[2]);

		return aliStrings;
	}

	private String[] mergeAlignments(String[] ali1, String[] ali2, int o) {
		String s1 = ali1[0].concat(ali2[0].substring(o));
		String s2 = ali1[1].concat(ali2[1].substring(o));
		String[] ali = { s1, s2 };
		return ali;
	}

	private String deriveSEQ(String[] ali) {
		String seq = ali[0].replaceAll("-", "");
		return seq;
	}

	private String deriveMD(String[] ali) {
		StringBuffer md = new StringBuffer("");
		char lastType = '-';
		int num = 0;
		for (int i = 0; i < ali[0].length(); i++) {
			char type = getMDType(ali[0].charAt(i), ali[1].charAt(i));
			if (type == 'M')
				num++;
			else if (type != 'I') {
				if (num != 0) {
					md = md.append(num);
					num = 0;
				}
				if (type == 'X') {
					if (lastType == 'D')
						md.append('0');
					md.append(ali[1].charAt(i));
				}
				if (type == 'D') {
					if (lastType != 'D')
						md.append('^');
					md.append(ali[1].charAt(i));
				}
			}
			if (type != 'I')
				lastType = type;
		}
		if (lastType == 'M')
			md = md.append(num);
		return md.toString();
	}

	private char getMDType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		if (c1 != c2)
			return 'X';
		return 'M';
	}

	private String deriveCIGAR(String[] ali) {
		StringBuffer cigar = new StringBuffer("");
		char lastType = getCIGARType(ali[0].charAt(0), ali[1].charAt(0));
		int num = 1;
		for (int i = 1; i < ali[0].length(); i++) {
			char type = getCIGARType(ali[0].charAt(i), ali[1].charAt(i));
			if (type == lastType)
				num++;
			else {
				cigar = cigar.append(num + "" + lastType);
				num = 1;
			}
			lastType = type;
		}
		cigar = cigar.append(num + "" + lastType);
		return cigar.toString();
	}

	private char getCIGARType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		return 'M';
	}

}
