package pipeline.init;

import io.Fasta_Writer;
import io.Fastq_Reader;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;

public class Shotgun {

	private HashMap<String, Integer> readToOrigLength = new HashMap<String, Integer>();

	public File run(File fastq, File out, int length, int step) {

		String name = fastq.getName().split("\\.")[0];
		File outFile = new File(out.getAbsolutePath() + File.separator + name + "_" + length + "_" + step + ".fna");
		outFile.delete();

		Vector<Read> long_reads = new Fastq_Reader().read(fastq);
		for (Read r : long_reads)
			readToOrigLength.put(r.getId(), r.getLength());

		System.out.println("Shredding " + long_reads.size() + " reads...");
		Fasta_Writer writer = new Fasta_Writer();
		int last_p = 0;
		for (Read r : long_reads) {
			String seq = r.getSeq();
			int counter = 1;
			for (int i = 0; i < seq.length() - length; i += step) {
				int j = i + length < seq.length() ? i + length : seq.length();
				String subseq = seq.substring(i, j);
				String id = r.getId() + ":" + (counter++);
				Read shortRead = new Read(id, subseq, "");
				writer.write(shortRead, outFile);
			}

			int p = (int) Math.round((double) (long_reads.indexOf(r) / (double) long_reads.size()) * 100.);
			if (p != last_p && p % 10 == 0) {
				System.out.println("OUTPUT>" + p + "% (" + long_reads.indexOf(r) + "/" + long_reads.size()
						+ ") of the reads shredded.");
				last_p = p;
			}

		}

		return outFile;

	}

	public HashMap<String, Integer> getReadToOrigLength() {
		return readToOrigLength;
	}

	public void setReadToOrigLength(HashMap<String, Integer> readToLength) {
		this.readToOrigLength = readToLength;
	}

}
