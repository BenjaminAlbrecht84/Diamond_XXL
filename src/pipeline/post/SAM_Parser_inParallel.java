package pipeline.post;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.SparseString;
import util.ScoringMatrix;

public class SAM_Parser_inParallel {

	private CountDownLatch latch;
	private ConcurrentHashMap<String, ReadHits> readMap;

	private AtomicInteger parsedLines = new AtomicInteger(0);
	private int last_p = 0, numOfLines;

	public ConcurrentHashMap<String, ReadHits> parse_Hits(File sam_file, ScoringMatrix matrix, HitRun_Rater scorer, int cores) {

		readMap = new ConcurrentHashMap<String, ReadHits>();

		try {

			// creating multiple sam_parserin respect to available cores
			numOfLines = countLines(sam_file);
			int chunk = (int) Math.ceil((double) numOfLines / (double) cores);
			chunk = chunk > 100000 ? 100000 : chunk;
			Vector<Parser> allParser = generateParser(sam_file, chunk, matrix);

			System.out.println("Parsing " + numOfLines + " lines...");

			// running sam_parsers in parallel
			latch = new CountDownLatch(allParser.size());
			ExecutorService executor = Executors.newFixedThreadPool(cores);
			for (Parser thread : allParser)
				executor.execute(thread);

			// awaiting termination
			try {
				latch.await();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			executor.shutdown();

			System.out.println("OUTPUT>" + 100 + "% (" + parsedLines + "/" + numOfLines + ") of the lines parsed.");

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return readMap;
	}

	private synchronized void reportProgress(int p) {
		p = ((int) Math.floor((double) p / 10.)) * 10;
		if (p != 100 && p != last_p) {
			System.out.println("OUTPUT>" + p + "% (" + parsedLines + "/" + numOfLines + ") of the lines parsed.");
			last_p = p;
		}
	}

	public Vector<Parser> generateParser(File file, int chunk, ScoringMatrix matrix) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(file));
		Vector<Parser> allParser = new Vector<Parser>();
		try {
			byte[] c = new byte[1024];
			int count = 0;
			int readChars = 0;
			long filePointer = 0;
			allParser.add(new Parser(file, matrix, filePointer, chunk));
			while ((readChars = is.read(c)) != -1) {
				for (int i = 0; i < readChars; ++i) {
					filePointer++;
					if (c[i] == '\n') {
						count++;
						if (count % (chunk + 1) == 0)
							allParser.add(new Parser(file, matrix, filePointer, chunk));
					}
				}
				break;
			}
			return allParser;
		} finally {
			is.close();
		}
	}

	public int countLines(File file) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(file));
		try {
			byte[] c = new byte[1024];
			int count = 0;
			int readChars = 0;
			boolean empty = true;
			while ((readChars = is.read(c)) != -1) {
				empty = false;
				for (int i = 0; i < readChars; ++i) {
					if (c[i] == '\n') {
						++count;
					}
				}
			}
			return (count == 0 && !empty) ? 1 : count;
		} finally {
			is.close();
		}
	}

	public class Parser implements Runnable {

		private File sam_file;
		private ScoringMatrix matrix;
		private int chunk;
		private long start;

		public Parser(File sam_file, ScoringMatrix matrix, long start, int chunk) {
			this.sam_file = sam_file;
			this.matrix = matrix;
			this.start = start;
			this.chunk = chunk;
		}

		public void run() {

			try {

				HashMap<String, ReadHits> localReadMap = new HashMap<String, ReadHits>();

				RandomAccessFile raf = new RandomAccessFile(sam_file, "r");
				raf.seek(start);
				long file_pointer = raf.getFilePointer();

				int readCounter = 0;
				String l;
				while ((l = raf.readLine()) != null) {

					if (localReadMap.keySet().size() == 1)
						break;

					if (readCounter > chunk)
						break;

					if (!l.startsWith("@") && !l.isEmpty()) {

						readCounter++;

						// START parsing row
						// **************************************

						// String[] columns = l.split("\t");
						String[] columns = mySplit(l, '\t');

						// shredding index
						String[] id_split = mySplit(columns[0], ':');
						int read_num = Integer.parseInt(id_split[id_split.length - 1]);

						// read_id without shredding index
						String read_id = columns[0].substring(0, columns[0].indexOf(":" + read_num));

						// starting position in the reference
						int ref_start = Integer.parseInt(columns[3]);

						// aligned query sequence without '*' at the end
						String seq = columns[9];
						seq = seq.charAt(seq.length() - 1) == '*' ? seq.substring(0, seq.length() - 1) : seq;

						// cigar sequence
						String cigar = columns[5];

						// everything else
						int bitScore = -1, rawScore = -1, frame = -1, query_start = -1, ref_length = -1;
						for (String c : columns) {
							// bit score
							if (c.startsWith("AS:i:"))
								bitScore = Integer.parseInt(mySplit(c, ':')[2]);
							// raw score
							else if (c.startsWith("ZR:i:"))
								rawScore = Integer.parseInt(mySplit(c, ':')[2]);
							// frame
							else if (c.startsWith("ZF:i:"))
								frame = Integer.parseInt(mySplit(c, ':')[2]);
							// starting position in the query
							else if (c.startsWith("ZS:i:"))
								query_start = Integer.parseInt(mySplit(c, ':')[2]);
							// length of the whole reference sequence
							else if (c.startsWith("ZL:i:"))
								ref_length = Integer.parseInt(mySplit(c, ':')[2]);
						}

						// derive refEnd from CIGAR
						int ref_end = ref_start - 1;
						Matcher matcher = Pattern.compile("[0-9]+[MD]+").matcher(cigar);
						while (matcher.find()) {
							String match = matcher.group();
							ref_end += Integer.parseInt(match.substring(0, match.length() - 1));
						}

						// derive queryLength from CIGAR
						int query_length = 0;
						matcher = Pattern.compile("[0-9]+[MI]+").matcher(cigar);
						while (matcher.find()) {
							String match = matcher.group();
							query_length += Integer.parseInt(match.substring(0, match.length() - 1));
						}

						// gi of the matched subsequence
						SparseString gi = new SparseString(columns[2]);
						// int gi =
						// Integer.parseInt(columns[2].split("\\|")[1]);

						// END parsing row
						// **************************************

						// Object[] scoring_result = new
						// ReconstructAlignment(matrix).run(seq, cigar,
						// mdString);
						// int[] alignmentScores = (int[]) scoring_result[0];
						// BitSet query_insertions = (BitSet) scoring_result[1];
						// BitSet query_deletions = (BitSet) scoring_result[2];
						// int aliScoresSum = (int) scoring_result[3];

						read_num = frame < 0 ? -read_num : read_num;

						System.out.println(read_num + " " + ref_start + " " + ref_end + " " + bitScore + " " + rawScore + " " + file_pointer + " "
								+ query_start + " " + ref_length + " " + query_length);

						// Storing Hit
						Hit h = new Hit(read_num, ref_start, ref_end, bitScore, rawScore, file_pointer, null, query_start, ref_length, query_length,
								-1);
						if (!localReadMap.containsKey(read_id))
							localReadMap.put(read_id, new ReadHits());
						localReadMap.get(read_id).add(h, gi, frame);

						if (readCounter % 1000 == 0) {
							int p = (int) Math.round(((double) parsedLines.addAndGet(1000) / (double) numOfLines) * 100.);
							reportProgress(p);
						}

					}

					// file_pointer = raf.getFilePointer();
					file_pointer += l.length() + 1;

				}

				raf.close();
				addHits(localReadMap);

			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			latch.countDown();

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

	private synchronized void addHits(HashMap<String, ReadHits> localReadMap) {
		for (String read_id : localReadMap.keySet()) {
			if (!readMap.containsKey(read_id))
				readMap.put(read_id, new ReadHits());
			ReadHits hits = localReadMap.get(read_id);
			for (SparseString gi : hits.getHitMap().keySet()) {
				for (int frame : hits.getHitMap().get(gi).keySet()) {
					for (Hit h : hits.getHitMap().get(gi).get(frame))
						readMap.get(read_id).add(h, gi, frame);
				}
			}
		}
	}

}
