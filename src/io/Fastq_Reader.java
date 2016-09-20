package io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import pipeline.init.Read;

public class Fastq_Reader {

	public Vector<Read> read(File fastq_File) {
		Vector<Read> reads = new Vector<Read>();
		try {
			BufferedReader buf = new BufferedReader(new FileReader(fastq_File));
			String line, id = "";
			boolean readSequence = false;
			StringBuffer seq = new StringBuffer("");
			StringBuffer qual = new StringBuffer("");
			while ((line = buf.readLine()) != null) {
				if (line.startsWith("@") || line.startsWith(">")) {
					if (seq.length() != 0 && !id.isEmpty())
						reads.add(new Read(id, seq.toString(), qual.toString()));
					seq = new StringBuffer("");
					qual = new StringBuffer("");
					id = line.substring(1).split(" ")[0];
					readSequence = true;
				} else if (line.startsWith("+")) {
					readSequence = false;
				} else if (readSequence) {
					seq = seq.append(line);
				} else {
					qual = qual.append(line);
				}
			}
			if (seq.length() != 0 && !id.isEmpty())
				reads.add(new Read(id, seq.toString(), qual.toString()));
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//		Vector<Read> filteredReads = new Vector<Read>();
//		for (int i = 366; i < 367; i++) {
//			System.out.println(reads.get(i).getId());
//			System.out.println(reads.get(i).getSeq());
//			filteredReads.add(reads.get(i));
//		}
//		return filteredReads;

		 return reads;
		 
	}

	public ConcurrentHashMap<String, Long> parseReadIDs(File fastq_File) {
		ConcurrentHashMap<String, Long> readToPointer = new ConcurrentHashMap<String, Long>();
		try {
			RandomAccessFile raf = new RandomAccessFile(fastq_File, "r");
			ByteBuffer buffer = ByteBuffer.allocate(1024);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			int readChars = 0;
			boolean parseID = false;
			StringBuffer buf = new StringBuffer();
			long pointer = 0;
			int lineCounter = 0;
			while ((readChars = raf.read(buffer.array())) != -1) {
				for (int i = 0; i < readChars; i++) {
					char c = (char) buffer.get(i);
					if (parseID && (c == ' ' || c == '\n')) {
						parseID = false;
					} else if (parseID) {
						buf = buf.append(c);
					} else if ((c == '@' && lineCounter == 0) || c == '>') {
						parseID = true;
					}
					if (c == '\n' && buf.length() != 0) {
						readToPointer.put(buf.toString(), pointer + i);
						buf = new StringBuffer();
					}
					if (c == '\n')
						lineCounter = lineCounter == 3 ? 0 : ++lineCounter;
				}
				pointer = raf.getFilePointer();
			}
			raf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return readToPointer;
	}
}
