package io;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.concurrent.ConcurrentLinkedQueue;

import io.daa.DAA_Hit;
import io.daa.DAA_Reader;
import pipeline.post.Hit;

public class SAM_Writer {

	private File sam_file;
	private DAA_Reader daaReader;

	public SAM_Writer(File sam_file, DAA_Reader daaReader) {
		this.sam_file = sam_file;
		this.daaReader = daaReader;
	}

	public synchronized void run(ConcurrentLinkedQueue<Hit> hits) {

		try {

			RandomAccessFile raf = new RandomAccessFile(sam_file, "rw");
			RandomAccessFile rafSAM = new RandomAccessFile(sam_file, "r");
			RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;

			try {

				long filePointer = sam_file.length();
				raf.seek(filePointer);

				for (Hit h : hits) {

					String[] aliStrings = h.getAlignmentStrings(null);

					String read_id, gi;
					if (daaReader != null) {
						DAA_Hit hit = daaReader.parseHit(rafDAA, h.getFile_pointer(), h.getAccessPoint());
						read_id = hit.getQueryName();
						gi = hit.getReferenceName();
					} else {
						rafSAM.seek(h.getFile_pointer());
						String[] columns = rafSAM.readLine().split("\t");
						read_id = columns[0];
						gi = columns[2];
					}
					
					StringBuffer buf = new StringBuffer("");
					buf = buf.append(read_id + "\t");
					buf = buf.append(0 + "\t");
					buf = buf.append(gi + "\t");
					buf = buf.append(h.getRef_start() + "\t");
					buf = buf.append(255 + "\t");
					buf = buf.append(aliStrings[0] + "\t");
					buf = buf.append("*\t");
					buf = buf.append("0\t");
					buf = buf.append("0\t");
					buf = buf.append(aliStrings[1] + "\t");
					buf = buf.append("*\t");
					buf = buf.append("AS:i:" + h.getBitScore() + "\t");
					buf = buf.append("ZL:i:" + h.getRef_length() + "\t");
					buf = buf.append("ZR:i:" + h.getRawScore() + "\t");
					buf = buf.append("ZS:i:" + h.getQuery_start() + "\t");
					buf = buf.append("MD:Z:" + aliStrings[2] + "\t");
					buf = buf.append("\n");

					// freeing memory
					h.setMetaInfo(null);
					h.freeMemory();

					String out = buf.toString();
					raf.write(out.getBytes());

					h.setFilePointer(filePointer);
					h.setAccessPoint(null);
					filePointer += out.length();

				}

			} finally {
				raf.close();
				rafSAM.close();
				if (rafDAA != null)
					rafDAA.close();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
