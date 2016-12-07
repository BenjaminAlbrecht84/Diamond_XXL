package io.debug;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import pipeline.post.Hit_Run;

public class RejectedWriter {

	private File rej_file;
	private boolean append = false;

	public RejectedWriter(File rej_file) {
		this.rej_file = rej_file;
	}

	public synchronized void run(Vector<Hit_Run> hitRuns) {

		try {

			BufferedWriter buf = new BufferedWriter(new FileWriter(rej_file, append));
			try {
				for (Hit_Run run : hitRuns) {
					buf.write(run.getReadID() + "\t");
					buf.write(run.getGi() + "\t");
					buf.write(run.getHitRun().size() + "\t");
					buf.write(run.getSumScore() + "\t");
					buf.write(run.getEValue() + "\t");
					buf.write(run.getCoverge() + "\t");
					buf.write("\n");
				}
			} finally {
				buf.close();
			}

			append = true;

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
