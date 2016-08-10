package io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Vector;

import pipeline.init.Read;
import util.Mutator;

public class Fastq_ErrorReader {

	public Vector<Read> read(File fastq_File, Mutator mut, double error_rate) {
		Vector<Read> reads = new Vector<Read>();
		try {
			BufferedReader buf = new BufferedReader(new FileReader(fastq_File));
			String line, id = "";
			boolean readSequence = false;
			StringBuffer seq = new StringBuffer("");
			StringBuffer qual = new StringBuffer("");
			while ((line = buf.readLine()) != null) {
				if (line.startsWith("@")) {
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
				if (seq.length() != 0 && !id.isEmpty())
					reads.add(new Read(id, mut.insert_insertions(seq.toString(), error_rate), qual.toString()));
			buf.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return reads;
	}

}
