package pipeline.post.mode_one;

import io.Dmnd_IndexReader;
import io.Fastq_Reader;
import io.daa.DAA_Reader;

import java.io.File;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.HitRun_Writer;
import pipeline.post.HitToSamConverter;
import pipeline.post.Hit_Run;
import pipeline.post.Hit.HitType;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.ReconstructAlignment;
import util.ScoringMatrix;
import util.frameshiftAligner.Frameshift_Alignment;
import util.frameshiftAligner.Frameshift_Alignment.AliMode;

public class Alignment_Completer {

	private CountDownLatch latch;
	private ConcurrentHashMap<String, Long> readToPointer;
	private Dmnd_IndexReader dmndReader;
	private File queryFile, refFile, samFile;
	private DAA_Reader daaReader;
	private ConcurrentHashMap<Integer, Character> indexToAA;
	private ScoringMatrix scoringMatrix;
	private double lambda, k;
	private HitRun_Rater hitRunRater;
	private int step;
	private HitToSamConverter samConverter;
	private HitRun_Writer runWriter;
	private Vector<Hit_Run> allRuns;
	private double maxEValue;
	private int minSumScore;

	private Vector<Alignment_Thread> aliThreads;
	private int last_p, totalNumberOfRuns;
	private AtomicInteger completedRuns;

	public void run(Vector<Hit_Run> allRuns, File queryFile, File refFile, File samFile, DAA_Reader daaReader, ScoringMatrix scoringMatrix,
			double lambda, double k, int cores, HitRun_Rater hitRunRater, int step, HitToSamConverter samConverter, HitRun_Writer runWriter,
			double maxEValue, int minSumScore) {

		this.allRuns = allRuns;
		this.queryFile = queryFile;
		this.refFile = refFile;
		this.samFile = samFile;
		this.daaReader = daaReader;
		this.scoringMatrix = scoringMatrix;
		this.k = k;
		this.lambda = lambda;
		this.hitRunRater = hitRunRater;
		this.step = step;
		this.samConverter = samConverter;
		this.runWriter = runWriter;
		this.maxEValue = maxEValue;
		this.minSumScore = minSumScore;

		completedRuns = new AtomicInteger(0);
		totalNumberOfRuns = allRuns.size();

		// initializing aa mapper
		String aaString = "ARNDCQEGHILKMFPSTWYVBJZX";
		indexToAA = new ConcurrentHashMap<Integer, Character>();
		for (int i = 0; i < aaString.length(); i++)
			indexToAA.put(i, aaString.charAt(i));

		// pre-computing file pointers
		readToPointer = new Fastq_Reader().parseReadIDs(queryFile);
		dmndReader = new Dmnd_IndexReader(refFile);
		dmndReader.createIndex();

		System.out.println("Closing " + allRuns.size() + " alignments...");

		// sorting all hit runs decreasingly after its sum score and reference
		// coverage
		Collections.sort(allRuns, new HitRun_Comparator());

		// alternately distributing runs on different alignment threads
		Vector<Vector<Hit_Run>> runSubsets = new Vector<Vector<Hit_Run>>();
		for (int i = 0; i < cores; i++)
			runSubsets.add(new Vector<Hit_Run>());
		aliThreads = new Vector<Alignment_Thread>();
		for (int i = 0; i < allRuns.size(); i++)
			runSubsets.get(i % cores).add(0, allRuns.get(i));
		for (int i = 0; i < cores; i++)
			aliThreads.add(new Alignment_Thread(runSubsets.get(i)));

		// running scorer in parallel
		latch = new CountDownLatch(aliThreads.size());
		ExecutorService executor = Executors.newFixedThreadPool(cores);
		for (Alignment_Thread thread : aliThreads)
			executor.execute(thread);

		// awaiting termination
		try {
			latch.await();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		executor.shutdown();

		int p = (int) Math.round(((double) completedRuns.get() / (double) totalNumberOfRuns) * 100.);
		System.out.println("OUTPUT>Finally, " + p + "% (" + completedRuns.get() + "/" + totalNumberOfRuns + ") of the alignments closed.");

	}

	private void reportCompletedRuns(Vector<Hit_Run> runs) {

		// writing all hits
		if (samConverter != null) {
			Vector<Hit> allHits = new Vector<Hit>();
			for (Hit_Run run : runs)
				allHits.addAll(run.getHitRun());
			samConverter.run(allHits);
		}

		// writing runs
		runWriter.run(runs);

	}

	private void reportProgess(Vector<Hit_Run> runs) {
		int p = (int) Math.round(((double) completedRuns.addAndGet(runs.size()) / (double) totalNumberOfRuns) * 100.);
		p = ((int) Math.floor((double) p / 1.)) * 1;
		if (p != 100 && p != last_p && p % 1 == 0) {
			System.out.println("OUTPUT>" + p + "% (" + completedRuns.get() + "/" + totalNumberOfRuns + ") of the most promising alignments closed.");
			last_p = p;
		}
	}

	public class Alignment_Thread implements Runnable {

		private Vector<Hit_Run> runs;

		public Alignment_Thread(Vector<Hit_Run> runs) {
			this.runs = runs;
		}

		@Override
		public void run() {

			try {

				RandomAccessFile rafSAM = new RandomAccessFile(samFile, "r");
				RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;

				try {

					int chunkSize = 10000;
					Vector<Hit_Run> subset = new Vector<Hit_Run>();

					while (!runs.isEmpty()) {

						for (int i = runs.size() - 1; i >= 0; i--) {
							Hit_Run run = runs.elementAt(i);
							subset.add(run);
							if (subset.size() % chunkSize == 0)
								break;
						}

						// loading GI Sequences
						Vector<Integer> gIs = new Vector<Integer>();
						for (Hit_Run run : subset) {
							int gi = run.getGi();
							if (!gIs.contains(gi))
								gIs.add(gi);
						}
						HashMap<Integer, String> giToSeq = getAASequences(dmndReader.getGILocations(gIs));

						// loading ReadID Sequences
						Vector<String> readIDs = new Vector<String>();
						for (Hit_Run run : subset) {
							String readID = run.getReadID();
							if (!readIDs.contains(readID))
								readIDs.add(readID);
						}
						HashMap<String, String> readIDToSeq = loadReadIDSequences(readIDs);

						// closing alignment gaps
						Vector<Hit_Run> highScoringRuns = new Vector<Hit_Run>();
						for (int j = 0; j < subset.size(); j++) {
							Hit_Run run = subset.get(j);
							closeAliGaps(run, giToSeq, readIDToSeq, rafSAM, rafDAA);
							if (run.getEValue() < maxEValue && run.getSumScore() > minSumScore)
								highScoringRuns.add(run);
						}

						// reporting processed runs
						reportCompletedRuns(highScoringRuns);
						reportProgess(subset);

						// freeing memory
						for (Hit_Run run : subset) {
							for (Hit h : run.getHitRun()) {
								h.freeMemory();
								h.setMetaInfo(null);
								h = null;
							}
							run.freeMemory();
							allRuns.remove(run);
						}
						highScoringRuns.clear();
						giToSeq.clear();
						readIDToSeq.clear();
						readIDs.clear();
						gIs.clear();

						// ready for processing next subset
						for (int i = 0; i < subset.size(); i++)
							runs.remove(runs.size() - 1);
						subset = new Vector<Hit_Run>();

					}

				} finally {
					rafSAM.close();
					if (rafDAA != null)
						rafDAA.close();
				}

			} catch (Exception e) {
				e.printStackTrace();
			}

			latch.countDown();

		}

		private HashMap<Integer, String> getAASequences(HashMap<Integer, Long> giToPointer) {

			HashMap<Integer, String> giToSeq = new HashMap<Integer, String>();
			Vector<Integer> sortedGIs = new Vector<Integer>();
			sortedGIs.addAll(giToPointer.keySet());
			Collections.sort(sortedGIs);

			try {
				RandomAccessFile raf = new RandomAccessFile(refFile, "r");
				try {
					for (int gi : sortedGIs) {
						raf.seek(giToPointer.get(gi));
						ByteBuffer buffer = ByteBuffer.allocate(1024);
						buffer.order(ByteOrder.LITTLE_ENDIAN);
						int readChars = 0;
						boolean doBreak = false;
						StringBuffer aaSeq = new StringBuffer();
						StringBuffer aaInt = new StringBuffer();
						while ((readChars = raf.read(buffer.array())) != -1 && !doBreak) {
							for (int r = 1; r < readChars; r++) {
								int aaIndex = (int) buffer.get(r);
								if (doBreak = (aaIndex == -1))
									break;
								aaSeq = aaSeq.append(indexToAA.get(aaIndex));
								aaInt = aaInt.append(aaIndex + ",");
							}
						}
						giToSeq.put(gi, aaSeq.toString());
					}
				} finally {
					raf.close();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			return giToSeq;
		}

		private void closeAliGaps(Hit_Run run, HashMap<Integer, String> giToSeq, HashMap<String, String> readIDToSeq, RandomAccessFile rafSAM,
				RandomAccessFile rafDAA) {

			// System.out.println("new Run: " + run.getReadID() + " " +
			// run.getGi() + " " + run.getFrameDirection());

			Frameshift_Alignment aligner = new Frameshift_Alignment(scoringMatrix, 20);
			Vector<Vector<Hit>> closingHits = new Vector<Vector<Hit>>();

			// loading query and reference sequences
			String query = readIDToSeq.get(run.getReadID());
			if (run.getFrameDirection() == Frame_Direction.Negativ)
				query = reverseComplementString(query);
			String ref = giToSeq.get(run.getGi());
			Frame_Direction frameDir = run.getFrameDirection();

			// closing first gap *********************************

			Hit h = run.getHitRun().get(0);
			int refGapStart = 0;
			int refGapEnd = h.getRef_start();

			int offset = Math.abs(h.getId() * step);
			int queryGapStart = 0;
			int queryGapEnd = frameDir == Frame_Direction.Positiv ? (h.getQuery_start() + offset) : query.length() - h.getQuery_start() - offset;

			int queryGapLength = queryGapEnd - queryGapStart + 1;
			int refGapLength = (refGapEnd - refGapStart + 1) * 3;
			queryGapStart = 2 * refGapLength < queryGapLength ? queryGapEnd - refGapLength : queryGapStart;

			Vector<Hit> closingFrameHits = null;
			if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {
				String subQuery = reverse(query.substring(queryGapStart, queryGapEnd));
				String subRef = reverse(ref.substring(refGapStart, refGapEnd));
				String[] aliResultRev = aligner.run(subQuery, subRef, AliMode.FREESHIFT);
				String[] aliResult = { reverse(aliResultRev[0]), reverse(aliResultRev[1]), reverse(aliResultRev[2]) };
				closingFrameHits = generateFrameHits(aliResult, queryGapStart, queryGapEnd, refGapStart, refGapEnd, h, run.getFrameDirection(),
						query.length());
			}
			closingHits.add(closingFrameHits);

			// closing middle gap *********************************

			queryGapStart = queryGapEnd + (h.getQuery_length() * 3);
			refGapStart = h.getRef_end();
			for (int i = 1; i < run.getHitRun().size(); i++) {
				closingFrameHits = null;
				if (h.getRef_end() < run.getHitRun().get(i).getRef_start()) {
					queryGapStart = queryGapEnd + (h.getQuery_length() * 3);
					refGapStart = h.getRef_end();
					h = run.getHitRun().get(i);
					refGapEnd = h.getRef_start();
					offset = Math.abs(h.getId() * step);
					queryGapEnd = frameDir == Frame_Direction.Positiv ? (h.getQuery_start() + offset) : query.length() - h.getQuery_start() - offset;
					closingFrameHits = null;
					if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {
						String subQuery = query.substring(queryGapStart, queryGapEnd);
						String subRef = ref.substring(refGapStart, refGapEnd);
						String[] aliResult = aligner.run(subQuery, subRef, AliMode.GLOBAL);
						closingFrameHits = generateFrameHits(aliResult, queryGapStart, queryGapEnd, refGapStart, refGapEnd, h,
								run.getFrameDirection(), query.length());
					}
				} else if (run.getHitRun().get(i).getRef_end() > h.getRef_end()) {
					h = run.getHitRun().get(i);
					offset = Math.abs(h.getId() * step);
					queryGapEnd = frameDir == Frame_Direction.Positiv ? (h.getQuery_start() + offset) : query.length() - h.getQuery_start() - offset;
					refGapStart = h.getRef_end();
				}
				closingHits.add(closingFrameHits);
			}

			// closing last gap *********************************

			queryGapStart = queryGapEnd + (h.getQuery_length() * 3);
			refGapStart = h.getRef_end();
			refGapEnd = ref.length();
			queryGapEnd = query.length();

			queryGapLength = queryGapEnd - queryGapStart + 1;
			refGapLength = (refGapEnd - refGapStart + 1) * 3;

			queryGapEnd = 2 * refGapLength < queryGapLength ? queryGapStart + refGapLength : queryGapEnd;

			closingFrameHits = null;
			if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {
				String subQuery = query.substring(queryGapStart, queryGapEnd);
				String subRef = giToSeq.get(run.getGi()).substring(refGapStart, refGapEnd);
				String[] aliResult = aligner.run(subQuery, subRef, AliMode.FREESHIFT);
				closingFrameHits = generateFrameHits(aliResult, queryGapStart, queryGapEnd, refGapStart, refGapEnd, h, run.getFrameDirection(),
						query.length());
			}
			closingHits.add(closingFrameHits);

			// ***************************************************

			// for (Hit hBef : run.getHitRun())
			// hBef.print("Before: ");

			// inserting closingHits into run
			int pos = closingHits.size() - 1;
			for (int i = run.getHitRun().size(); i >= 0; i--) {
				if (closingHits.get(pos) != null) {
					for (int j = closingHits.get(pos).size() - 1; j >= 0; j--) {
						Hit closingHit = closingHits.get(pos).get(j);
						run.getHitRun().insertElementAt(closingHit, i);
					}
				}
				pos--;
			}

			// for (Hit hBef : run.getHitRun())
			// hBef.print("After: ");

			// updating hitRun
			run.update(hitRunRater, rafSAM, rafDAA);
			run.setCompleted(true);

		}

		private String reverse(String s) {
			return new StringBuffer(s).reverse().toString();
		}

		private String reverseComplementString(String s) {
			String rev = new StringBuffer(s).reverse().toString();
			String revComp = rev.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
			return revComp;
		}

		private Vector<Hit> generateFrameHits(String[] aliResult, int queryStart, int queryEnd, int refStart, int refEnd, Hit hit,
				Frame_Direction frame_Direction, int totalQueryLength) {

			Vector<Hit> closingHits = new Vector<Hit>();

			StringBuffer[] subAliResult = new StringBuffer[3];
			subAliResult[0] = new StringBuffer();
			subAliResult[1] = new StringBuffer();
			int qStart = queryStart, qLength = 0, rStart = refStart, rLength = 0;

			String frameIDs = aliResult[2];
			Integer lastFrameID = 1;
			for (int i = 0; i <= frameIDs.length(); i++) {

				Integer frameID = i < frameIDs.length() ? Character.getNumericValue(frameIDs.charAt(i)) : null;

				if (i == 0)
					qStart = qStart + frameID - 1;

				if (i < frameIDs.length() && (i == 0 || lastFrameID == frameID)) {

					// extending frame hit
					char q = aliResult[0].charAt(i);
					char r = aliResult[1].charAt(i);
					subAliResult[0] = subAliResult[0].append(q);
					subAliResult[1] = subAliResult[1].append(r);
					qLength = q != '-' ? qLength + 3 : qLength;
					rLength = r != '-' ? rLength + 1 : rLength;

				} else {

					// reporting frame hit
					String queryAli = subAliResult[0].toString();
					String refAli = subAliResult[1].toString();
					String[] subAli = { queryAli, refAli };
					Hit h = generateHit(subAli, qStart, qStart + qLength, rStart, rStart + rLength, hit, lastFrameID, frame_Direction,
							totalQueryLength);
					closingHits.add(h);

					if (i < frameIDs.length()) {

						// computing offset at query start
						int offset = frameID - lastFrameID;

						// resetting parameters for recording next framehit
						qStart = qStart + qLength + offset;
						rStart = rStart + rLength;
						qLength = 0;
						rLength = 0;
						subAliResult[0] = new StringBuffer();
						subAliResult[1] = new StringBuffer();
						i--;

					}

				}
				lastFrameID = frameID;
			}

			return closingHits;

		}

		private Hit generateHit(String[] aliResult, int queryStart, int queryEnd, int refStart, int refEnd, Hit hit, int frame,
				Frame_Direction frame_Direction, int totalQueryLength) {

			// initializing closing hit
			int score = 0;
			for (int s : scoringMatrix.cmpAlignmentScores(aliResult[0], aliResult[1]))
				score += s;
			int rawScore = score;
			int bitScore = (int) Math.round((lambda * (double) rawScore - Math.log(k)) / Math.log(2));
			int refLength = hit.getRef_length();
			int queryLength = (queryEnd - queryStart + 1) / 3;
			queryStart = frame_Direction == Frame_Direction.Positiv ? queryStart : totalQueryLength - queryStart;
			Hit h = new Hit(-1, refStart, refEnd, bitScore, rawScore, hit.getFile_pointer(), hit.getAccessPoint(), queryStart, refLength,
					queryLength);
			h.setHitType(HitType.Synthetic);

			// adding metaInfo
			frame = frame_Direction == Frame_Direction.Positiv ? frame : -frame;
			Object[] metaInfo = { queryStart, frame };
			h.setMetaInfo(metaInfo);

			// deriving alignment strings
			String[] ali = { aliResult[0], aliResult[1] };
			String[] aliStrings = new String[3];
			aliStrings[0] = deriveCIGAR(ali);
			aliStrings[1] = deriveSEQ(ali);
			aliStrings[2] = deriveMD(ali);
			h.copyAliStrings(aliStrings);

			return h;
		}

		private HashMap<String, String> loadReadIDSequences(Vector<String> readIDs) {
			HashMap<String, String> readIDToSeq = new HashMap<String, String>();
			try {
				RandomAccessFile raf = new RandomAccessFile(queryFile, "r");
				for (String readID : readIDs) {
					Long filePointer = readToPointer.get(readID);
					raf.seek(filePointer + 1);

					ByteBuffer buffer = ByteBuffer.allocate(1024);
					buffer.order(ByteOrder.LITTLE_ENDIAN);
					int readChars = 0;
					StringBuffer buf = new StringBuffer();
					boolean doBreak = false;
					while (((readChars = raf.read(buffer.array())) != -1) && !doBreak) {
						for (int i = 0; i < readChars; i++) {
							char c = (char) buffer.get(i);
							if (c == '\n') {
								doBreak = true;
								break;
							}
							buf = buf.append(c);
						}
					}
					String seq = buf.toString();
					readIDToSeq.put(readID, seq);
				}
				raf.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			return readIDToSeq;
		}

		// private HashMap<String, String> loadReadIDSequences(Vector<String>
		// readIDs) {
		// HashMap<String, String> readIDToSeq = new HashMap<String, String>();
		// try {
		// RandomAccessFile raf = new RandomAccessFile(queryFile, "r");
		// for (String readID : readIDs) {
		// Long filePointer = readToPointer.get(readID);
		// raf.seek(filePointer + 1);
		// String seq = raf.readLine();
		// readIDToSeq.put(readID, seq);
		// }
		// raf.close();
		// } catch (Exception e) {
		// e.printStackTrace();
		// }
		// return readIDToSeq;
		// }

		private String deriveSEQ(String[] ali) {
			String seq = ali[0].replaceAll("-", "");
			return seq;
		}

		private String deriveMD(String[] ali) {
			StringBuffer md = new StringBuffer("");
			char lastType = '-';
			int num = 0;
			for (int i = 0; i < ali[0].length(); i++) {
				char type = getMDType(ali[0].charAt(i), ali[1].charAt(i));
				if (type == 'M')
					num++;
				else if (type != 'I') {
					if (num != 0) {
						md = md.append(num);
						num = 0;
					}
					if (type == 'X') {
						if (lastType == 'D')
							md.append('0');
						md.append(ali[1].charAt(i));
					}
					if (type == 'D') {
						if (lastType != 'D')
							md.append('^');
						md.append(ali[1].charAt(i));
					}
				}
				if (type != 'I')
					lastType = type;
			}
			if (lastType == 'M')
				md = md.append(num);
			return md.toString();
		}

		private char getMDType(char c1, char c2) {
			if (c1 == '-')
				return 'D';
			if (c2 == '-')
				return 'I';
			if (c1 != c2)
				return 'X';
			return 'M';
		}

		private String deriveCIGAR(String[] ali) {
			StringBuffer cigar = new StringBuffer("");
			char lastType = getCIGARType(ali[0].charAt(0), ali[1].charAt(0));
			int num = 1;
			for (int i = 1; i < ali[0].length(); i++) {
				char type = getCIGARType(ali[0].charAt(i), ali[1].charAt(i));
				if (type == lastType)
					num++;
				else {
					cigar = cigar.append(num + "" + lastType);
					num = 1;
				}
				lastType = type;
			}
			cigar = cigar.append(num + "" + lastType);
			return cigar.toString();
		}

		private char getCIGARType(char c1, char c2) {
			if (c1 == '-')
				return 'D';
			if (c2 == '-')
				return 'I';
			return 'M';
		}

	}

	public class HitRun_Comparator implements Comparator<Hit_Run> {

		// @Override
		// public int compare(Hit_Run r1, Hit_Run r2) {
		// if (r1.getGi() < r2.getGi())
		// return -1;
		// if (r1.getGi() > r2.getGi())
		// return 1;
		// return 0;
		// }

		@Override
		public int compare(Hit_Run r1, Hit_Run r2) {
			// int s1 = (int) Math.round(Math.log(Math.pow((double)
			// r1.getScore(), (double) r1.getCoverge())));
			// int s2 = (int) Math.round(Math.log(Math.pow((double)
			// r2.getScore(), (double) r2.getCoverge())));
			int s1 = r1.getSumScore() * r1.getCoverge();
			int s2 = r2.getSumScore() * r2.getCoverge();
			// int s1 = r1.getCoverge();
			// int s2 = r2.getCoverge();
			if (s1 > s2)
				return -1;
			if (s1 < s2)
				return 1;
			return 0;
		}

	}

}
