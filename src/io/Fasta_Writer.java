package io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import pipeline.init.Read;

public class Fasta_Writer {

	public void write(Vector<Read> reads, File out) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(out, true));
			for (Read r : reads) {
				writer.write(">" + r.getId());
				writer.newLine();
				writer.write(r.getSeq());
				writer.newLine();
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
