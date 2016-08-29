package pipeline.post;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.math.BigInteger;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.GroupLayout.Alignment;

import io.daa.DAA_Hit;
import io.daa.DAA_Reader;
import util.CompressAlignment;
import util.ScoringMatrix;

public class HitToSamConverter {

	private File sam_file;
	private DAA_Reader daaReader;
	private int step;
	private BigInteger dbSize;
	private HashMap<String, Integer> readToOrigLength;
	private ScoringMatrix matrix;
	private File out;

	public HitToSamConverter(File sam_file, DAA_Reader daaReader, int step, BigInteger dbSize, HashMap<String, Integer> readToOrigLength,
			ScoringMatrix matrix, File out) {
		this.sam_file = sam_file;
		this.daaReader = daaReader;
		this.step = step;
		this.dbSize = dbSize;
		this.readToOrigLength = readToOrigLength;
		this.matrix = matrix;
		this.out = out;
		writeHeader();
	}

	private void writeHeader() {
		try {

			StringBuffer header = new StringBuffer("@PG\tPN:DIAMOND\n");
			header = header.append("@mm\tBlastX" + "\n");
			header = header.append("@CO\tBlastX-like alignments\n");
			header = header.append(
					"@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate\n");

			FileWriter writer = new FileWriter(out);
			writer.write(header.toString());
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void run(List<Hit> hits) {

		try {

			StringBuffer buf = new StringBuffer("");
			RandomAccessFile rafSAM = new RandomAccessFile(sam_file, "r");
			RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;

			try {

				for (Hit hit : hits) {

					if (hit.getAccessPoint() == null) {

						rafSAM.seek(hit.getFile_pointer());
						String l = rafSAM.readLine();
						String[] columns = l.split("\t");

						// parsing shredding index
						String[] id_split = columns[0].split(":");
						int read_num = Integer.parseInt(id_split[id_split.length - 1]);

						// parsing read_id without shredding index
						String read_id = columns[0].substring(0, columns[0].indexOf(":" + read_num));

						// parsing frame and query_start
						Integer frame = null, query_start = -1;
						if (hit.getMetaInfo() == null) {
							query_start = (read_num) * step;
							for (String c : columns) {
								// frame of the query
								if (c.startsWith("ZF:i:"))
									frame = Integer.parseInt(c.split(":")[2]);
								// starting position in the query
								else if (c.startsWith("ZS:i:"))
									query_start += Integer.parseInt(c.split(":")[2]);
							}
						} else {
							query_start = (Integer) hit.getMetaInfo()[0];
							frame = (Integer) hit.getMetaInfo()[1];
						}

						// parsing gi/ref
						String gi = columns[2];

						// determining ali strings
						String[] aliStrings = hit.getAlignmentStrings(rafSAM);

						// computing eValue
						BigInteger m = dbSize;
						double n = readToOrigLength.get(read_id);
						double sPrime = hit.getBitScore();
						double eValue = m.doubleValue() * n * Math.pow(2, -sPrime);

						// computing percent identity
						double matches = 0, aliLength = hit.getAlignmentScores(matrix, rafSAM).length;
						Matcher matcher = Pattern.compile("[0-9]+").matcher(aliStrings[2]);
						while (matcher.find())
							matches += Integer.parseInt(matcher.group());
						int identity = (int) Math.floor((matches / aliLength) * 100.);

						buf = buf.append(read_id + "\t");
						buf = buf.append(0 + "\t");
						buf = buf.append(gi + "\t");
						buf = buf.append(hit.getRef_start() + "\t");
						buf = buf.append(255 + "\t");
						buf = buf.append(aliStrings[0] + "\t");
						buf = buf.append("*\t");
						buf = buf.append("0\t");
						buf = buf.append("0\t");
						buf = buf.append(aliStrings[1]);
						buf = buf.append("*\t");
						buf = buf.append("AS:i:" + hit.getBitScore() + "\t");
						buf = buf.append("ZL:i:" + hit.getRef_length() + "\t");
						buf = buf.append("ZR:i:" + hit.getRawScore() + "\t");
						DecimalFormat dfEValue = new DecimalFormat("0.##E0");
						buf = buf.append("ZE:f:" + dfEValue.format(eValue) + "\t");
						buf = buf.append("ZI:i:" + identity + "\t");
						if (frame != null)
							buf = buf.append("ZF:i:" + frame + "\t");
						buf = buf.append("ZS:i:" + query_start + "\t");
						buf = buf.append("MD:Z:" + aliStrings[2] + "\t");
						buf = buf.append("\n");

					} else {

						DAA_Hit h = daaReader.parseHit(rafDAA, hit.getFile_pointer(), hit.getAccessPoint());

						String read_id = h.getQueryName();
						String gi = h.getReferenceName();

						// parsing read number
						String[] id_split = mySplit(read_id, ':');
						read_id = id_split[0];
						int read_num = Integer.parseInt(id_split[id_split.length - 1]);

						// parsing frame and query_start
						Integer frame = null, query_start = -1;
						if (hit.getMetaInfo() == null) {
							int offset = read_num * step;
							query_start = h.getQueryStart() + 1 + offset;
							frame = h.getFrame();
						} else {
							query_start = (Integer) hit.getMetaInfo()[0];
							frame = (Integer) hit.getMetaInfo()[1];
						}

						// determining ali strings
						String[] aliStrings = new CompressAlignment().run(h.getAlignment());

						// computing eValue
						BigInteger m = dbSize;
						double n = readToOrigLength.get(read_id);
						double sPrime = hit.getBitScore();
						double eValue = m.doubleValue() * n * Math.pow(2, -sPrime);

						// computing percent identity
						String[] alignment = h.getAlignment();
						double matches = 0, aliLength = alignment.length;
						for (int i = 0; i < alignment[0].length(); i++)
							matches = alignment[0].charAt(i) == alignment[1].charAt(i) ? matches + 1 : matches;
						int identity = (int) Math.floor((matches / aliLength) * 100.);

						buf = buf.append(read_id + "\t");
						buf = buf.append(0 + "\t");
						buf = buf.append(gi + "\t");
						buf = buf.append(hit.getRef_start() + "\t");
						buf = buf.append(255 + "\t");
						buf = buf.append(aliStrings[0] + "\t");
						buf = buf.append("*\t");
						buf = buf.append("0\t");
						buf = buf.append("0\t");
						buf = buf.append(aliStrings[1]);
						buf = buf.append("*\t");
						buf = buf.append("AS:i:" + hit.getBitScore() + "\t");
						buf = buf.append("ZL:i:" + hit.getRef_length() + "\t");
						buf = buf.append("ZR:i:" + hit.getRawScore() + "\t");
						DecimalFormat dfEValue = new DecimalFormat("0.##E0");
						buf = buf.append("ZE:f:" + dfEValue.format(eValue) + "\t");
						buf = buf.append("ZI:i:" + identity + "\t");
						buf = buf.append("ZF:i:" + frame + "\t");
						buf = buf.append("ZS:i:" + query_start + "\t");
						buf = buf.append("MD:Z:" + aliStrings[2] + "\t");
						buf = buf.append("\n");

					}

				}

			} finally {
				rafSAM.close();
				if (rafDAA != null)
					rafDAA.close();
			}

			writeInFile(buf);

		} catch (Exception e) {
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

	private String[] mySplit(String s, char c) {
		List<String> words = new ArrayList<String>();
		int pos = 0, end;
		while ((end = s.indexOf(c, pos)) >= 0) {
			words.add(s.substring(pos, end));
			pos = end + 1;
		}
		if (pos < s.length())
			words.add(s.substring(pos, s.length()));
		String[] entries = words.toArray(new String[words.size()]);
		return entries;

	}

}
