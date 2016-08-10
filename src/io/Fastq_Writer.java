package io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import pipeline.init.Read;


public class Fastq_Writer {

	public void write(Vector<Read> reads, File out) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(out));
			for (Read r : reads) {
				writer.write("@" + r.getId());
				writer.newLine();
				writer.write(r.getSeq());
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(r.getQual());
				writer.newLine();
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void write(Read r, File out, boolean doAppend) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(out, doAppend));
			writer.write("@" + r.getId());
			writer.newLine();
			writer.write(r.getSeq());
			writer.newLine();
			writer.write("+");
			writer.newLine();
			writer.write(r.getQual());
			writer.newLine();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
