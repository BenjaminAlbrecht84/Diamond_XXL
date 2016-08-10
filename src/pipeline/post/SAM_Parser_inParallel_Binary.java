package pipeline.post;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
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

import util.ScoringMatrix;

public class SAM_Parser_inParallel_Binary {

	private CountDownLatch latch;
	private ConcurrentHashMap<String, ReadHits> readMap;

	private AtomicInteger allParsedLines = new AtomicInteger(0);
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

			System.out.println("OUTPUT>" + 100 + "% (" + allParsedLines + "/" + numOfLines + ") of the lines parsed.");

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return readMap;
	}

	private synchronized void reportProgress(int p) {
		p = ((int) Math.floor((double) p / 10.)) * 10;
		if (p != 100 && p != last_p && p % 1 == 0) {
			System.out.println("OUTPUT>" + p + "% (" + allParsedLines + "/" + numOfLines + ") of the lines parsed.");
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
				// break;
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
				long file_pointer = start;
				long pointer = start;

				int parsedLines = 0;

				byte[] buffer = new byte[1024 * 1024];
				boolean readingHeader = true;
				int readChars = 0, colNumber = 0;
				char lastChar = '\n';

				StringBuffer buf = new StringBuffer();
				int gi = -1, read_num = -1, frame = -1, ref_start = -1, bitScore = -1, rawScore = -1, query_start = -1, ref_length = -1;
				String cigar = "", read_id = "", seq = "", mdz = "";

				StringBuffer allEntries = new StringBuffer("");
				while ((readChars = raf.read(buffer)) != -1) {

					if (parsedLines > chunk)
						break;

					for (int i = 0; i < readChars; i++) {

						char c = (char) buffer[i];

						if (c == '\n')
							parsedLines++;

						if (parsedLines > chunk)
							break;

						if (readingHeader && lastChar == '\n' && c != '@') {
							readingHeader = false;
							pointer += i;
						}

						if (!readingHeader) {

							switch (c) {
							case '\n':

								read_num = frame < 0 ? -read_num : read_num;

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

								// generating hit
								Hit h = new Hit(read_num, ref_start, ref_end, bitScore, rawScore, pointer, null, query_start, ref_length,
										query_length);

								// storing compressed alignment strings (needs
								// too much memory!!!)
								// String[] aliStrings = new String[3];
								// aliStrings[0] = cigar;
								// aliStrings[1] = seq;
								// aliStrings[2] = mdz;
								// h.setAliStrings(aliStrings);

								// storing hit
								if (!localReadMap.containsKey(read_id))
									localReadMap.put(read_id, new ReadHits());
								localReadMap.get(read_id).add(h, gi, frame);

								// reporting progress
								if (parsedLines % 1000 == 0) {
									int p = (int) Math.round(((double) allParsedLines.addAndGet(1000) / (double) numOfLines) * 100.);
									reportProgress(p);
								}

								buf = new StringBuffer();
								pointer = file_pointer + i + 1;
								colNumber = 0;

								break;
							case '\t':
								String entry = buf.toString();
								allEntries = allEntries.append(entry + "\t");
								switch (colNumber) {
								case (0):
									String[] id_split = mySplit(entry, ':');
									read_id = id_split[0];
									read_num = Integer.parseInt(id_split[id_split.length - 1]);
									break;
								case (2):
									gi = Integer.parseInt(mySplit(entry, '|')[1]);
									break;
								case (3):
									ref_start = Integer.parseInt(entry);
									break;
								case (5):
									cigar = entry;
									break;
								case (9):
									seq = entry;
									seq = seq.charAt(seq.length() - 1) == '*' ? seq.substring(0, seq.length() - 1) : seq;
									break;
								case (11):
									bitScore = Integer.parseInt(mySplit(entry, ':')[2]);
									break;
								case (13):
									ref_length = Integer.parseInt(mySplit(entry, ':')[2]);
									break;
								case (14):
									rawScore = Integer.parseInt(mySplit(entry, ':')[2]);
									break;
								case (17):
									frame = Integer.parseInt(mySplit(entry, ':')[2]);
									break;
								case (18):
									query_start = Integer.parseInt(mySplit(entry, ':')[2]);
									break;
								case (19):
									mdz = mySplit(entry, ':')[2];
									break;
								}
								colNumber++;
								buf = new StringBuffer();
								break;
							default:
								buf = buf.append(c);
							}

						}

						lastChar = c;

					}

					file_pointer = raf.getFilePointer();

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

	private synchronized void addHit(String read_id, Hit h, int gi, int frame) {
		if (!readMap.containsKey(read_id))
			readMap.put(read_id, new ReadHits());
		readMap.get(read_id).add(h, gi, frame);
	}

	private synchronized void addHits(HashMap<String, ReadHits> localReadMap) {
		for (String read_id : localReadMap.keySet()) {
			if (!readMap.containsKey(read_id))
				readMap.put(read_id, new ReadHits());
			ReadHits hits = localReadMap.get(read_id);
			for (int gi : hits.getHitMap().keySet()) {
				for (int frame : hits.getHitMap().get(gi).keySet()) {
					for (Hit h : hits.getHitMap().get(gi).get(frame))
						readMap.get(read_id).add(h, gi, frame);
				}
			}
		}
	}

}
