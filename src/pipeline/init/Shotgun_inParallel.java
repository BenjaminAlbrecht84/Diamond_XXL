package pipeline.init;

import io.Fasta_Writer;
import io.Fastq_Reader;

import java.io.File;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class Shotgun_inParallel {

	private HashMap<String, Integer> readToOrigLength = new HashMap<String, Integer>();
	private CountDownLatch latch;

	private AtomicInteger shred_counter = new AtomicInteger(0);
	private int last_p = 0, totalReads;

	public File run(File fastq, File out, int length, int step, int cores) {
		
		String name = fastq.getName().split("\\.")[0];
		File outFile = new File(out.getAbsolutePath() + File.separator + name + "_" + length + "_" + step + ".fna");
		outFile.delete();

		Vector<Read> long_reads = new Fastq_Reader().read(fastq);
		for (Read r : long_reads)
			readToOrigLength.put(r.getId(), r.getLength());

		totalReads = long_reads.size();
		
		System.out.println("Shredding " + totalReads + " read...");

		Fasta_Writer writer = new Fasta_Writer();
		int chunkSize = (int) Math.ceil((double) long_reads.size() / (double) cores);
		Vector<Read> read_chunk = new Vector<Read>();
		Vector<Shredder> allShredder = new Vector<Shredder>();
		for (int i = 0; i < long_reads.size(); i++) {
			read_chunk.add(long_reads.get(i));
			if ((i != 0 && i % chunkSize == 0) || long_reads.size() == 1) {
				allShredder.add(new Shredder(read_chunk, length, step, writer, outFile));
				read_chunk = new Vector<Read>();
			}
		}
		allShredder.add(new Shredder(read_chunk, length, step, writer, outFile));

		latch = new CountDownLatch(allShredder.size());
		ExecutorService executor = Executors.newFixedThreadPool(cores);
		for (Shredder thread : allShredder)
			executor.execute(thread);

		try {
			latch.await();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("OUTPUT>100% (" + shred_counter + "/" + totalReads + ") of the reads shredded.");

		executor.shutdown();

		return outFile;

	}

	public class Shredder implements Runnable {

		private Vector<Read> long_reads = new Vector<Read>();
		private int length, step;
		private Fasta_Writer writer;
		private File outFile;

		public Shredder(Vector<Read> long_reads, int length, int step, Fasta_Writer writer, File outFile) {
			this.long_reads = long_reads;
			this.length = length;
			this.step = step;
			this.writer = writer;
			this.outFile = outFile;
		}

		public void run() {
			for (Read r : long_reads) {
				String seq = r.getSeq();
				int counter = 0;
				Vector<Read> readBuffer = new Vector<Read>();
				for (int i = 0; i < seq.length(); i += step) {
					
					int j = i + length < seq.length() ? i + length : seq.length();
					String subseq = seq.substring(i, j);
					String id = r.getId() + ":" + (counter++);
					readBuffer.add(new Read(id, subseq, ""));
					
					if(readBuffer.size() % 1000 == 0){
						writeReads(writer, readBuffer, outFile);
						readBuffer = new Vector<Read>();
					}
					
					if (j == seq.length())
						break;
				}
				writeReads(writer, readBuffer, outFile);
				int p = (int) Math.round(((double) shred_counter.incrementAndGet() / (double) totalReads) * 100.);
				reportProgress(p);
			}
			latch.countDown();
		}

	}

	private synchronized void writeReads(Fasta_Writer writer, Vector<Read> shortReads, File outFile) {
		for (Read r : shortReads)
			writer.write(r, outFile);
	}

	private synchronized void reportProgress(int p) {
		if (p != 100 && p != last_p && p % 10 == 0) {
			System.out.println("OUTPUT>" + p + "% (" + shred_counter + "/" + totalReads + ") of the reads shredded.");
			last_p = p;
		}
	}

	public HashMap<String, Integer> getReadToOrigLength() {
		return readToOrigLength;
	}

	public void setReadToOrigLength(HashMap<String, Integer> readToLength) {
		this.readToOrigLength = readToLength;
	}

}
