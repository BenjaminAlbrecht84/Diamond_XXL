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
import util.ScoringMatrix;

public class HitRun_Writer {

	private File out, samFile;
	private DAA_Reader daaReader;

	public HitRun_Writer(File out, File samFile, DAA_Reader daaReader) {
		this.out = out;
		this.samFile = samFile;
		this.daaReader = daaReader;
		writeHeader();
	}

	private void writeHeader() {

		StringBuffer buf = new StringBuffer("");
		buf = buf.append("@PG\t PN:DIAMOND_XXL\n");
		buf = buf.append("@mm\t BlastX\n");
		buf = buf.append("@CO\t BlastX-like alignments\n");
		buf = buf.append("@CO\t AS: sumScore, ZS: rawScore, ZE: expected, ZC: coverage, ZQ: query coverlength, ZL: reference length\n");

		try {
			FileWriter fW = new FileWriter(out, false);
			fW.write(buf.toString());
			fW.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public synchronized void run(Vector<Hit_Run> runs) {

		try {

			StringBuffer buf = new StringBuffer("");
			RandomAccessFile rafSAM = new RandomAccessFile(samFile, "r");
			RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;
			Collections.sort(runs, new HitRunComparator());
			for (Hit_Run run : runs) {
				Hit h = run.getHitRun().firstElement();
				String reference = h.getFile_pointer() != -1 ? getReference(h, h.getFile_pointer(), rafSAM, rafDAA) : "gi|" + run.getGi();
				DecimalFormat dfEValue = new DecimalFormat("0.##E0");
				buf = buf.append(run.getReadID() + "\t" + reference + "\t AS:i:" + run.getSumScore() + "\t ZR:i:" + run.getRawScore() + "\t ZE:f:"
						+ dfEValue.format(run.getEValue()) + "\t ZC:i:" + run.getCoverge() + "\t ZQ:i:" + run.getRunLength() + "\t ZL:i:"
						+ run.getHitRun().get(0).getRef_length() + "\n");
			}
			rafSAM.close();
			if (rafDAA != null)
				rafDAA.close();

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
