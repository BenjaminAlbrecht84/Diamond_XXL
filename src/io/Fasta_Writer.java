package io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import pipeline.init.Read;

public class Fasta_Writer {

	public void write(Read r, File out) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(out, true));
			writer.write(">" + r.getId());
			writer.newLine();
			writer.write(r.getSeq());
			writer.newLine();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
