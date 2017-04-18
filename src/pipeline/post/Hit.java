package pipeline.post;

import java.io.RandomAccessFile;
import java.util.BitSet;

import io.daa.DAA_Hit;
import io.daa.DAA_Reader;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.AlignmentEvaluater;
import util.CompressAlignment;
import util.ReconstructAlignment;
import util.ScoringMatrix;

public class Hit {

	public enum HitType {
		Real, Synthetic, Merged
	}

	private HitType hitType = HitType.Real;
	private int frame = Integer.MIN_VALUE;

	private int id, ref_start, ref_end, bitScore, rawScore, query_start, query_length, ref_length;
	private Long file_pointer;
	private Integer accessPoint;
	private int subjectID = -1;

	private int[] alignmentScores;
	private BitSet query_insertions;
	private BitSet query_deletions;

	public int[] aliStats;
	public String[] aliStrings;
	private Object[] metaInfo = null;

	public Hit(int id, int ref_start, int ref_end, int bitScore, int rawScore, long file_pointer, Integer accessPoint, int query_start,
			int ref_length, int query_length, int[] alignmentScores, BitSet query_insertions, BitSet query_deletions, int subjectID) {

		this.id = new Integer(id);
		this.ref_start = new Integer(ref_start);
		this.ref_end = new Integer(ref_end);
		this.bitScore = new Integer(bitScore);
		this.rawScore = new Integer(rawScore);
		this.file_pointer = new Long(file_pointer);
		this.accessPoint = accessPoint != null ? new Integer(accessPoint) : null;
		this.query_start = new Integer(query_start);
		this.ref_length = new Integer(ref_length);
		this.query_length = new Integer(query_length);
		this.subjectID = subjectID;

		this.alignmentScores = new int[alignmentScores.length];
		for (int i = 0; i < alignmentScores.length; i++)
			this.alignmentScores[i] = new Integer(alignmentScores[i]);
		this.query_insertions = (BitSet) query_insertions.clone();
		this.query_deletions = (BitSet) query_deletions.clone();

	}

	public Hit(int id, int ref_start, int ref_end, int bitScore, int rawScore, long file_pointer, Integer accessPoint, int query_start,
			int ref_length, int query_length, int subjectID) {
		this.id = new Integer(id);
		this.ref_start = new Integer(ref_start);
		this.ref_end = new Integer(ref_end);
		this.bitScore = new Integer(bitScore);
		this.rawScore = new Integer(rawScore);
		this.file_pointer = new Long(file_pointer);
		this.accessPoint = accessPoint != null ? new Integer(accessPoint) : null;
		this.query_start = new Integer(query_start);
		this.ref_length = new Integer(ref_length);
		this.query_length = new Integer(query_length);
		this.subjectID = subjectID;
	}

	public Hit(Hit h) {
		this.id = new Integer(h.getId());
		this.ref_start = new Integer(h.getRef_start());
		this.ref_end = new Integer(h.getRef_end());
		this.bitScore = new Integer(h.getBitScore());
		this.rawScore = new Integer(h.getRawScore());
		this.file_pointer = new Long(h.getFile_pointer());
		this.accessPoint = h.getAccessPoint() != null ? new Integer(h.getAccessPoint()) : null;
		this.query_start = new Integer(h.getQuery_start());
		this.ref_length = new Integer(h.getRef_length());
		this.query_length = new Integer(h.getQuery_length());
		this.subjectID = new Integer(h.getSubjectID());
		this.frame = new Integer(h.getFrame());
	}

	public int numOfQueryInsertions(int l, int r) {
		if (l <= r) {
			int i = l, counter = 0;
			while (i - counter < r) {
				if (query_insertions.get(i))
					counter++;
				i++;
			}
			return counter;
		} else {
			int i = l, counter = 0;
			while (i + counter > r) {
				if (query_insertions.get(i))
					counter++;
				i--;
			}
			return counter;
		}
	}

	public boolean isQueryInsertion(int pos) {
		return query_insertions.get(pos);
	}

	public int numOfQueryDeletions(int l, int r) {
		int i = l, counter = 0;
		while (i - counter < r) {
			if (query_deletions.get(i))
				counter++;
			i++;
		}
		return counter;
	}

	public int numOfQueryDeletionsFixed(int l, int r) {
		int i = l, counter = 0;
		while (i < r) {
			if (query_deletions.get(i))
				counter++;
			i++;
		}
		return counter;
	}

	public boolean isQueryDeletion(int pos) {
		return query_deletions.get(pos);
	}

	public int[] getAlignmentStats(ScoringMatrix matrix, RandomAccessFile raf) {
		if (aliStats != null)
			return aliStats;
		String[] aliStrings = getAlignmentStrings(raf);
		Object[] scoringResult = new ReconstructAlignment(matrix).run(aliStrings[1], aliStrings[0], aliStrings[2]);
		String[] alignment = { (String) scoringResult[4], (String) scoringResult[5] };
		aliStats = new AlignmentEvaluater().run(alignment, matrix);
		return aliStats;
	}

	public int[] getAlignmentStats(ScoringMatrix matrix, RandomAccessFile raf, DAA_Reader daaReader) {
		if (aliStats != null)
			return aliStats;
		DAA_Hit daaHit = daaReader.parseHit(raf, file_pointer, accessPoint);
		aliStats = new AlignmentEvaluater().run(daaHit.getAlignment(), matrix);
		return aliStats;
	}

	public int[] getAlignmentScores(ScoringMatrix matrix, RandomAccessFile raf) {
		if (alignmentScores != null)
			return alignmentScores;
		String[] aliStrings = getAlignmentStrings(raf);
		Object[] scoring_result = new ReconstructAlignment(matrix).run(aliStrings[1], aliStrings[0], aliStrings[2]);
		alignmentScores = (int[]) scoring_result[0];
		query_insertions = (BitSet) scoring_result[1];
		query_deletions = (BitSet) scoring_result[2];
		return alignmentScores;
	}

	public int[] getAlignmentScores(ScoringMatrix matrix, RandomAccessFile raf, DAA_Reader daaReader) {
		if (alignmentScores != null)
			return alignmentScores;
		String[] aliStrings = getAlignmentStrings(raf, daaReader);
		Object[] scoring_result = new ReconstructAlignment(matrix).run(aliStrings[1], aliStrings[0], aliStrings[2]);
		alignmentScores = (int[]) scoring_result[0];
		query_insertions = (BitSet) scoring_result[1];
		query_deletions = (BitSet) scoring_result[2];
		return alignmentScores;
	}

	public String[] getAlignmentStrings(RandomAccessFile raf, DAA_Reader daaReader) {
		if (getAliStrings() != null)
			return getAliStrings();
		try {
			DAA_Hit daaHit = daaReader.parseHit(raf, file_pointer, accessPoint);
			String[] aliStrings = new CompressAlignment().run(daaHit.getAlignment());
			return aliStrings;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public String[] getAlignmentStrings(RandomAccessFile raf) {
		if (getAliStrings() != null)
			return getAliStrings();
		try {
			raf.seek(file_pointer);
			String line = raf.readLine();
			String[] columns = line.split("\t");
			String[] aliStrings = new String[3];
			aliStrings[0] = columns[5];
			aliStrings[1] = columns[9];
			for (String c : columns) {
				if (c.startsWith("MD:Z:"))
					aliStrings[2] = c.split(":")[2];
			}
			return aliStrings;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public void setId(int id) {
		this.id = id;
	}

	public int getId() {
		return id;
	}

	public int getRef_start() {
		return ref_start;
	}

	public int getRef_end() {
		return ref_end;
	}

	public int getBitScore() {
		return bitScore;
	}

	public int getRawScore() {
		return rawScore;
	}

	public Long getFile_pointer() {
		return file_pointer;
	}

	public Integer getAccessPoint() {
		return accessPoint;
	}

	public int getQuery_start() {
		return query_start;
	}

	public void setQuery_start(int qStart) {
		this.query_start = qStart;
	}

	public int getQuery_length() {
		return query_length;
	}

	public int getRef_length() {
		return ref_length;
	}

	public int getAliLength() {
		return alignmentScores.length;
	}

	public void print(String prefix) {
		int queryEnd = frame > 0 ? (query_start + 3 * query_length) : (query_start - 3 * query_length);
		System.out.println(prefix + " " + id + "\tQB:[" + query_start + "," + queryEnd + "]\tRB:[" + ref_start + "," + ref_end
				+ "]\tRS: " + rawScore + "\tBS: " + bitScore + "\tFR: " + frame);
	}

	public String toString() {
		return (id + "\tQB:[" + query_start + "," + query_length + "]\tRB:[" + ref_start + "," + ref_end + "]\tRS: " + rawScore + "\tBS: " + bitScore
				+ "\tFR: " + frame);
	}

	public String[] getAliStrings() {
		return aliStrings;
	}

	public void copyAliStrings(String[] aliStrings) {
		this.aliStrings = new String[3];
		this.aliStrings[0] = new String(aliStrings[0]);
		this.aliStrings[1] = new String(aliStrings[1]);
		this.aliStrings[2] = new String(aliStrings[2]);
	}

	public void setMetaInfo(Object[] metaInfo) {
		this.metaInfo = metaInfo;
	}

	public Object[] getMetaInfo() {
		return metaInfo;
	}

	public void setHitType(HitType hitType) {
		this.hitType = hitType;
	}

	public HitType getHitType() {
		return hitType;
	}

	public void freeMemory() {
		aliStrings = null;
		alignmentScores = null;
		file_pointer = null;
	}

	public void setFilePointer(long file_pointer) {
		this.file_pointer = new Long(file_pointer);
	}

	public void setAccessPoint(Integer accessPoint) {
		this.accessPoint = accessPoint;
	}

	public int getFrame() {
		return frame;
	}

	public void setFrame(int frame) {
		this.frame = frame;
	}

	public int getSubjectID() {
		return subjectID;
	}

	public void setAlignmentStats(int[] aliStats) {
		this.aliStats = aliStats;
	}

	public BitSet getQueryInsertions() {
		return query_insertions;
	}

}
