package pipeline.post;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

import io.daa.DAA_Hit;
import io.daa.DAA_Reader;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.ScoringMatrix;

public class HitRun_Writer {

	private File out, samFile;
	private DAA_Reader daaReader;
	private ScoringMatrix matrix;

	public HitRun_Writer(File out, File samFile, DAA_Reader daaReader, ScoringMatrix matrix) {
		this.out = out;
		this.samFile = samFile;
		this.daaReader = daaReader;
		this.matrix = matrix;
		writeHeader();
	}

	private void writeHeader() {

		StringBuffer buf = new StringBuffer("");
		buf = buf.append("@PG\t PN:DIAMOND_XXL\n");
		buf = buf.append("@mm\t BlastX\n");
		buf = buf.append("@CO\t BlastX-like alignments\n");
		buf = buf.append(
				"@CO\t qseqid sseqid length matches positves mismatches gapopens gaps qstart qend sstart send scoverage frames bitScore eValue\n");

		try {
			FileWriter fW = new FileWriter(out, false);
			fW.write(buf.toString());
			fW.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void run(Vector<Hit_Run> runs) {

		try {

			StringBuffer buf = new StringBuffer("");
			RandomAccessFile rafSAM = new RandomAccessFile(samFile, "r");
			RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;
//			Collections.sort(runs, new HitRunComparator());

			for (Hit_Run run : runs) {

				Hit hFirst = run.getHitRun().firstElement();
				Hit hLast = run.getHitRun().lastElement();
				String reference = hFirst.getFile_pointer() != -1 ? getReference(hFirst, hFirst.getFile_pointer(), rafSAM, rafDAA)
						: "gi|" + run.getGi();
				int[] aliStats = run.getAliStats();
				DecimalFormat dfEValue = new DecimalFormat("0.##E0");

				buf = buf.append(run.getReadID() + "\t");
				buf = buf.append(reference + "\t");
				buf = buf.append(aliStats[0] + "\t");
				buf = buf.append(aliStats[1] + "\t");
				buf = buf.append(aliStats[2] + "\t");
				buf = buf.append(aliStats[3] + "\t");
				buf = buf.append(aliStats[4] + "\t");
				buf = buf.append(aliStats[5] + "\t");
				buf = buf.append(hFirst.getQuery_start() + "\t");
				if (run.getFrameDirection() == Frame_Direction.Positiv)
					buf = buf.append((hLast.getQuery_start() + hLast.getQuery_length() * 3 - 1) + "\t");
				else
					buf = buf.append((hLast.getQuery_start() - hLast.getQuery_length() * 3 + 1) + "\t");
				buf = buf.append(hFirst.getRef_start() + "\t");
				buf = buf.append(hLast.getRef_end() + "\t");
				buf = buf.append(run.getCoverge() + "\t");
				buf = buf.append(run.getHitRun().size() + "\t");
				buf = buf.append(run.getSumScore() + "\t");
				buf = buf.append(dfEValue.format(run.getEValue()) + "\t");
				buf = buf.append("\n");

			}
			rafSAM.close();
			if (rafDAA != null)
				rafDAA.close();

			writeInFile(buf);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private synchronized void writeInFile(StringBuffer buf) {
		try {
			FileWriter fW = new FileWriter(out, true);
			fW.write(buf.toString());
			fW.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private String getReference(Hit h, long filePointer, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		try {
			if (h.getAccessPoint() != null) {
				DAA_Hit hit = daaReader.parseHit(rafDAA, h.getFile_pointer(), h.getAccessPoint());
				return hit.getReferenceName();
			}
			rafSAM.seek(h.getFile_pointer());
			String line = rafSAM.readLine();
			return line.split("\\s+")[2];
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "";
	}

	public class HitRunComparator implements Comparator<Hit_Run> {

		@Override
		public int compare(Hit_Run r1, Hit_Run r2) {
			int s1 = r1.getSumScore();
			int s2 = r2.getSumScore();
			if (s1 > s2)
				return -1;
			if (s1 < s2)
				return 1;
			return 0;
		}

	}

}
