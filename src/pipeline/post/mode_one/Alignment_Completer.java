package pipeline.post.mode_one;

import io.Dmnd_IndexReader;
import io.Fastq_Reader;
import io.daa.DAA_Reader;
import io.daa.DAA_Writer;
import io.debug.RejectedWriter;

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

import pipeline.post.Alignment_Merger;
import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.HitRun_Writer;
import pipeline.post.HitToSamConverter;
import pipeline.post.Hit_Run;
import pipeline.post.Hit.HitType;
import pipeline.post.HitRun_Linearizer;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import testArea.simulation.ReportAlignmentFiles;
import util.AA_Alphabet;
import util.AlignmentEvaluater;
import util.CodonTranslator;
import util.CompressAlignment;
import util.DAACompressAlignment;
import util.SparseString;
import util.ReconstructAlignment;
import util.ScoringMatrix;
import util.frameshiftAligner.banded.Banded_Frameshift_Alignment;
import util.frameshiftAligner.normal.Frameshift_Alignment;
import util.frameshiftAligner.sparse.Sparse_Frameshift_Alignment;

public class Alignment_Completer {

	static final double BLOCK_SIZE = 15000000;

	private CountDownLatch latch;
	private ConcurrentHashMap<String, Long> readToPointer;
	private Dmnd_IndexReader dmndReader;
	private File queryFile, refFile, samFile;
	private DAA_Reader daaReader;
	private ConcurrentHashMap<Integer, Character> indexToAA;
	private ScoringMatrix scoringMatrix;
	private double lambda, k;
	private HitRun_Rater hitRunRater;
	private int step, length;
	private HitToSamConverter samConverter;
	private DAA_Writer daaWriter;
	private HitRun_Writer runWriter;
	private Vector<Hit_Run> allRuns;
	private double maxEValue;
	private int minBitScore, minCoverage;
	private boolean useFilters, realign;

	private RejectedWriter rejectedWriter;
	private ReportAlignmentFiles aliReporter;

	private Vector<Alignment_Thread> aliThreads;
	private int last_p, totalNumberOfRuns;
	private AtomicInteger completedRuns, reportedRuns, rejectedRuns;

	public void run(Vector<Hit_Run> runs, File queryFile, File refFile, File samFile, DAA_Reader daaReader, ScoringMatrix scoringMatrix,
			double lambda, double k, int cores, HitRun_Rater hitRunRater, int step, int length, HitToSamConverter samConverter, DAA_Writer daaWriter,
			HitRun_Writer runWriter, double maxEValue, int minSumScore, int minCoverage, boolean useFilters, boolean realign, File rej_file_2,
			double blockSize, ReportAlignmentFiles aliReporter) {

		this.allRuns = runs;
		this.queryFile = queryFile;
		this.refFile = refFile;
		this.samFile = samFile;
		this.daaReader = daaReader;
		this.scoringMatrix = scoringMatrix;
		this.k = k;
		this.lambda = lambda;
		this.hitRunRater = hitRunRater;
		this.step = step;
		this.length = length;
		this.samConverter = samConverter;
		this.daaWriter = daaWriter;
		this.runWriter = runWriter;
		this.maxEValue = maxEValue;
		this.minBitScore = minSumScore;
		this.minCoverage = minCoverage;
		this.useFilters = useFilters;
		this.realign = realign;

		this.rejectedWriter = rej_file_2 != null ? new RejectedWriter(rej_file_2) : null;
		this.aliReporter = aliReporter;

		rejectedRuns = new AtomicInteger(0);
		reportedRuns = new AtomicInteger(0);
		completedRuns = new AtomicInteger(0);
		totalNumberOfRuns = allRuns.size();

		// initializing AA-Mapper
		String aaString = new AA_Alphabet().getAaString();
		indexToAA = new ConcurrentHashMap<Integer, Character>();
		for (int i = 0; i < aaString.length(); i++)
			indexToAA.put(i, aaString.charAt(i));

		// computing number of database chunks (30000000 sequences correspond to ~6GB)
		int chunkSize = (int) Math.round(BLOCK_SIZE * blockSize);
		int totalNumOfSeqs = new Dmnd_IndexReader(refFile).getNumberOfSequences();
		int totalIndexChunks = (int) Math.ceil((double) totalNumOfSeqs / (double) chunkSize);

		for (int indexChunk = 0; indexChunk < totalIndexChunks; indexChunk++) {

			// pre-computing file pointers
			readToPointer = new Fastq_Reader().parseReadIDs(queryFile);
			dmndReader = new Dmnd_IndexReader(refFile, indexChunk, totalIndexChunks);
			dmndReader.createIndex();

			long time = System.currentTimeMillis();
			System.out.println("STEP_6>Computing " + allRuns.size() + " alignment(s)...");

			// sorting all hit runs decreasingly after its sum score and reference coverage
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

			long runtime = (System.currentTimeMillis() - time) / 1000;
			int p = (int) Math.round(((double) completedRuns.get() / (double) totalNumberOfRuns) * 100.);
			System.out.println(
					"OUTPUT>" + p + "% (" + completedRuns.get() + "/" + totalNumberOfRuns + ") of the alignments computed. [" + runtime + "s]\n");

		}

		System.out.println("OUTPUT>" + reportedRuns.get() + " alignments reported!");

		if (daaWriter != null)
			daaWriter.finish();

	}

	private void reportCompletedRuns(Vector<Hit_Run> runs, HashMap<String, String> readIDToSeq) {

		// writing all hits in DAA file
		if (daaWriter != null)
			daaWriter.run(runs, readIDToSeq);

		// writing runs
		runWriter.run(runs);

		// updating reportedRuns counter
		reportedRuns.getAndAdd(runs.size());

	}

	private void reportProgess(int processedRuns) {
		int p = (int) Math.round(((double) completedRuns.addAndGet(processedRuns) / (double) totalNumberOfRuns) * 100.);
		p = ((int) Math.floor((double) p / 1.)) * 1;
		if (p != 100 && p != last_p && p % 1 == 0) {
			System.out.println(
					"OUTPUT>" + p + "% (" + completedRuns.get() + "/" + totalNumberOfRuns + ") of the most promising alignments computed...");
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

					int chunkSize = 5000;
					Vector<Hit_Run> subset = new Vector<Hit_Run>();

					while (!runs.isEmpty()) {

						for (int i = runs.size() - 1; i >= 0; i--) {
							Hit_Run run = runs.elementAt(i);
							subset.add(run);
							if (subset.size() % chunkSize == 0)
								break;
						}

						// loading GI Sequences
						Vector<SparseString> gIs = new Vector<SparseString>();
						for (Hit_Run run : subset) {
							SparseString gi = run.getGi();
							if (!gIs.contains(gi))
								gIs.add(gi);
						}
						HashMap<SparseString, String> giToSeq = getAASequences(gIs);

						// loading ReadID Sequences
						Vector<String> readIDs = new Vector<String>();
						for (Hit_Run run : subset) {
							String readID = run.getReadID();
							if (!readIDs.contains(readID))
								readIDs.add(readID);
						}
						HashMap<String, String> readIDToSeq = loadReadIDSequences(readIDs);

						// closing alignment gaps
						int processedRuns = 0;
						Vector<Hit_Run> highScoringRuns = new Vector<Hit_Run>();
						Vector<Hit_Run> lowScoringRuns = new Vector<Hit_Run>();
						for (int j = 0; j < subset.size(); j++) {
							Hit_Run run = subset.get(j);
							if (giToSeq.containsKey(run.getGi())) {

								// System.out.println("1-----------------------");
								// printRun(run, rafSAM, rafDAA);

								closeAliGaps(run, giToSeq, readIDToSeq, rafSAM, rafDAA, true);

								// System.out.println("2-----------------------");
								// printRun(run, rafSAM, rafDAA);
								// checkRun(run, rafSAM, rafDAA, readIDToSeq.get(run.getReadID()));

								mergeAlignments(run, rafSAM, rafDAA, readIDToSeq, giToSeq);

								// System.out.println("3-----------------------");
								// printRun(run, rafSAM, rafDAA);
								// checkRun(run, rafSAM, rafDAA, readIDToSeq.get(run.getReadID()));

								run.update(hitRunRater, rafSAM, rafDAA, daaReader, scoringMatrix);
								run.setCompleted(true);
								if (!useFilters || (run.getSumScore() > minBitScore && run.getCoverge() > minCoverage))
									highScoringRuns.add(run);
								else if (rejectedWriter != null)
									lowScoringRuns.add(run);
								processedRuns++;

								// System.out.println("4-----------------------");
								// printRun(run, rafSAM, rafDAA);
								// checkRun(run, rafSAM, rafDAA, readIDToSeq.get(run.getReadID()));

							}
						}

						// reporting processed runs
						reportCompletedRuns(highScoringRuns, readIDToSeq);
						reportProgess(processedRuns);
						if (aliReporter != null)
							aliReporter.run(highScoringRuns, rafSAM, rafDAA, readIDToSeq);

						// reporting rejected runs
						if (rejectedWriter != null)
							rejectedWriter.run(lowScoringRuns);
						rejectedRuns.addAndGet(lowScoringRuns.size());
						lowScoringRuns.clear();

						// freeing memory
						for (Hit_Run run : subset) {
							if (run.isCompleted()) {
								for (Hit h : run.getHitRun()) {
									h.freeMemory();
									h.setMetaInfo(null);
									h = null;
								}
								run.freeMemory();
								allRuns.remove(run);
							}
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

		private void mergeAlignments(Hit_Run run, RandomAccessFile rafSAM, RandomAccessFile rafDAA, HashMap<String, String> readIDToSeq,
				HashMap<SparseString, String> giToSeq) {
			// Hit h1 = run.getHitRun().get(0);
			boolean hitsMerged = false;
			for (int i = 1; i < run.getHitRun().size(); i++) {

				Hit h1 = run.getHitRun().get(i - 1);
				Hit h2 = run.getHitRun().get(i);
				int q1Start = h1.getQuery_start();
				int q1End = run.getFrameDirection() == Frame_Direction.Positiv ? q1Start + h1.getQuery_length() * 3
						: q1Start - h1.getQuery_length() * 3 - 1;
				int q2Start = h2.getQuery_start();

				if (h1.getFrame() == h2.getFrame() && h1.getRef_end() + 1 == h2.getRef_start() && (q1End + 1) == q2Start) {

					Hit h = new Alignment_Merger(hitRunRater, samFile, daaReader, step, length).mergeTwoHits(h1, h2, scoringMatrix, rafSAM, rafDAA,
							run, readIDToSeq, giToSeq, run.getFrameDirection());

					int pos = run.getHitRun().indexOf(h1);
					run.getHitRun().removeElement(h1);
					run.getHitRun().removeElement(h2);
					run.getHitRun().insertElementAt(h, pos);
					hitsMerged = true;
					break;
				}
			}
			if (hitsMerged)
				mergeAlignments(run, rafSAM, rafDAA, readIDToSeq, giToSeq);
		}

		private HashMap<SparseString, String> getAASequences(Vector<SparseString> gIs) {

			HashMap<SparseString, String> giToSeq = new HashMap<SparseString, String>();
			try {
				RandomAccessFile raf = new RandomAccessFile(refFile, "r");
				try {
					for (SparseString gi : gIs) {
						Long loc = dmndReader.getGILocation(gi);
						if (loc != null) {
							loc += 1;
							raf.seek(loc);
							ByteBuffer buffer = ByteBuffer.allocate(1024);
							buffer.order(ByteOrder.LITTLE_ENDIAN);
							int readChars = 0;
							boolean doBreak = false;
							StringBuffer aaSeq = new StringBuffer();
							while ((readChars = raf.read(buffer.array())) != -1 && !doBreak) {
								for (int r = 0; r < readChars; r++) {
									int aaIndex = (int) buffer.get(r);
									if (doBreak = (aaIndex == -1))
										break;
									if (indexToAA.containsKey(aaIndex))
										aaSeq = aaSeq.append(indexToAA.get(aaIndex));
									else
										aaSeq = aaSeq.append("X");
								}
							}
							giToSeq.put(gi, aaSeq.toString());
						}
					}
				} finally {
					raf.close();
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			return giToSeq;
		}

		private void closeAliGaps(Hit_Run run, HashMap<SparseString, String> giToSeq, HashMap<String, String> readIDToSeq, RandomAccessFile rafSAM,
				RandomAccessFile rafDAA, boolean firstRun) {

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
			int refGapEnd = h.getRef_start() - 1;

			int offset = h.getHitType() != HitType.Synthetic ? Math.abs(h.getId() * step) : 0;
			int queryGapStart = 0;
			int queryGapEnd = frameDir == Frame_Direction.Positiv ? (h.getQuery_start() - 1 + offset)
					: query.length() - h.getQuery_start() - 1 - offset;

			int queryGapLength = queryGapEnd - queryGapStart + 1;
			int refGapLength = (refGapEnd - refGapStart + 1) * 3;

			int queryGapStartOffset = 0;
			if (2 * refGapLength < queryGapLength) {
				queryGapStartOffset = queryGapEnd - refGapLength;
				queryGapStart = queryGapEnd - refGapLength;
			}
			// queryGapStart = 2 * refGapLength < queryGapLength ? queryGapEnd - refGapLength : queryGapStart;

			Vector<Hit> closingFrameHits = null;
			if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {

				String subQuery = query.substring(queryGapStart, queryGapEnd + 3);
				String subRef = ref.substring(refGapStart, refGapEnd);

				Object[] aliResult = new Banded_Frameshift_Alignment(scoringMatrix, scoringMatrix.getGapOpen()).run(subQuery, subRef,
						Banded_Frameshift_Alignment.AliMode.FREESHIFT_LEFT, 0.1, 0.1);
				// Object[] aliResult = new Frameshift_Alignment(scoringMatrix, 2 * scoringMatrix.getGapOpen()).run(subQuery, subRef,
				// Frameshift_Alignment.AliMode.FREESHIFT_LEFT);
				queryGapStart = ((int) (aliResult[5]) * 3) - (((int[]) aliResult[4])[0] * 3) + queryGapStartOffset;
				refGapStart = subRef.length() - ((int[]) aliResult[4])[1];

				closingFrameHits = generateFrameHits(aliResult, queryGapStart + 1, refGapStart + 1, h, run.getFrameDirection(), query.length(), run);
			}
			closingHits.add(closingFrameHits);

			// closing middle gap *********************************

			for (int i = 1; i < run.getHitRun().size(); i++) {
				closingFrameHits = null;
				if (h.getRef_end() < run.getHitRun().get(i).getRef_start()) {
					queryGapStart = queryGapEnd + (h.getQuery_length() * 3);
					refGapStart = h.getRef_end();
					h = run.getHitRun().get(i);
					refGapEnd = h.getRef_start() - 1;
					offset = h.getHitType() != HitType.Synthetic ? Math.abs(h.getId() * step) : 0;
					queryGapEnd = frameDir == Frame_Direction.Positiv ? (h.getQuery_start() + offset) : query.length() - h.getQuery_start() - offset;
					closingFrameHits = null;
					if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {
						String subQuery = query.substring(queryGapStart, queryGapEnd);
						String subRef = ref.substring(refGapStart, refGapEnd);
						Object[] aliResult = new Banded_Frameshift_Alignment(scoringMatrix, scoringMatrix.getGapOpen()).run(subQuery, subRef,
								Banded_Frameshift_Alignment.AliMode.GLOBAL, 0.1, 0.1);
						// Object[] aliResult = new Frameshift_Alignment(scoringMatrix, 2 * scoringMatrix.getGapOpen()).run(subQuery, subRef,
						// Frameshift_Alignment.AliMode.GLOBAL);
						closingFrameHits = generateFrameHits(aliResult, queryGapStart + 1, refGapStart + 1, h, run.getFrameDirection(),
								query.length(), run);
					}
				} else if (run.getHitRun().get(i).getRef_end() > h.getRef_end()) {
					h = run.getHitRun().get(i);
					offset = h.getHitType() != HitType.Synthetic ? Math.abs(h.getId() * step) : 0;
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
				Object[] aliResult = new Banded_Frameshift_Alignment(scoringMatrix, scoringMatrix.getGapOpen()).run(subQuery, subRef,
						Banded_Frameshift_Alignment.AliMode.FREESHIFT_RIGHT, 0.1, 0.1);
				// Object[] aliResult = new Frameshift_Alignment(scoringMatrix, 2 * scoringMatrix.getGapOpen()).run(subQuery, subRef,
				// Frameshift_Alignment.AliMode.FREESHIFT_RIGHT);
				closingFrameHits = generateFrameHits(aliResult, queryGapStart + 1, refGapStart + 1, h, run.getFrameDirection(), query.length(), run);
			}
			closingHits.add(closingFrameHits);

			// ***************************************************

			// inserting closingHits into run
			boolean synHitsAdded = false;
			int pos = closingHits.size() - 1;
			for (int i = run.getHitRun().size(); i >= 0; i--) {
				if (closingHits.get(pos) != null) {
					for (int j = closingHits.get(pos).size() - 1; j >= 0; j--) {
						Hit closingHit = closingHits.get(pos).get(j);
						run.getHitRun().insertElementAt(closingHit, i);
						synHitsAdded = true;
					}
				}
				pos--;
			}

			// updating hitRun
			Collections.sort(run.getHitRun(), new Hit_Comparator());

			// re-computing diamondHits **************************

			if (realign) {
				if (synHitsAdded) {

					// removing real hits
					Vector<Integer> toRemove = new Vector<Integer>();
					for (int i = 0; i < run.getHitRun().size(); i++) {
						if (run.getHitRun().get(i).getHitType() != HitType.Synthetic)
							toRemove.add(0, i);
					}

					// calling method recursively
					if (!toRemove.isEmpty()) {
						for (int i : toRemove)
							run.getHitRun().removeElementAt(i);
						closeAliGaps(run, giToSeq, readIDToSeq, rafSAM, rafDAA, false);
					}

				} else if (firstRun) {

					// calculating semi-global alignment
					h = run.getHitRun().get(0);
					refGapStart = 0;
					refGapEnd = ref.length();
					queryGapStart = 0;
					queryGapEnd = query.length();

					if (refGapStart < refGapEnd && queryGapStart < queryGapEnd - 4) {
						Object[] aliResult = new Banded_Frameshift_Alignment(scoringMatrix, scoringMatrix.getGapOpen()).run(query, ref,
								Banded_Frameshift_Alignment.AliMode.SEMI_GLOBAL, 1., 1.);
						queryGapStart = ((int[]) aliResult[4])[0] * 3;
						refGapStart = ((int[]) aliResult[4])[1];
						Vector<Hit> frameHits = generateFrameHits(aliResult, queryGapStart + 1, refGapStart + 1, h, run.getFrameDirection(),
								query.length(), run);
						run.getHitRun().clear();
						run.getHitRun().addAll(frameHits);
					}

				}
			} else {

				// adapting query start and frames
				for (Hit hit : run.getHitRun()) {
					if (hit.getHitType() != HitType.Synthetic) {
						int qStart = Math.abs(hit.getId()) * step + hit.getQuery_start() - 1;
						hit.setQuery_start(qStart);
						int frame = run.getFrameDirection() == Frame_Direction.Positiv ? (qStart % 3) + 1
								: -(((query.length() - (qStart + 1)) % 3) + 1);
						hit.setFrame(frame);
					}
				}

				// resolving overlaps
				new HitRun_Linearizer(daaReader, scoringMatrix, lambda, k).run(run, rafSAM, rafDAA, query);

			}

		}

		public class Hit_Comparator implements Comparator<Hit> {

			@Override
			public int compare(Hit h1, Hit h2) {
				if (h1.getRef_start() < h2.getRef_start())
					return -1;
				if (h1.getRef_start() > h2.getRef_start())
					return 1;
				return 0;
			}

		}

		private Vector<Hit> generateFrameHits(Object[] aliResult, int queryStart, int refStart, Hit hit, Frame_Direction frame_Direction,
				int totalQueryLength, Hit_Run run) {

			Vector<Hit> closingHits = new Vector<Hit>();

			StringBuffer[] subAliResult = new StringBuffer[3];
			subAliResult[0] = new StringBuffer();
			subAliResult[1] = new StringBuffer();
			int qStart = queryStart, qLength = 0, rStart = refStart, rLength = 0;

			String frameIDs = (String) aliResult[2];
			Integer lastFrameID = 1;
			for (int i = 0; i <= frameIDs.length(); i++) {

				Integer frameID = i < frameIDs.length() ? Character.getNumericValue(frameIDs.charAt(i)) : null;

				if (i == 0)
					qStart = qStart + frameID - 1;

				if (i < frameIDs.length() && (i == 0 || lastFrameID == frameID)) {

					// extending frame hit
					char q = ((String) aliResult[0]).charAt(i);
					char r = ((String) aliResult[1]).charAt(i);
					subAliResult[0] = subAliResult[0].append(q);
					subAliResult[1] = subAliResult[1].append(r);
					qLength = q != '-' ? qLength + 3 : qLength;
					rLength = r != '-' ? rLength + 1 : rLength;

				} else {

					// reporting frame hit
					String queryAli = subAliResult[0].toString();
					String refAli = subAliResult[1].toString();
					String[] subAli = { queryAli, refAli };

					if (rLength > 1) {
						Hit h = generateHit(subAli, qStart, qStart + qLength, rStart, rStart + rLength - 1, hit, lastFrameID, frame_Direction,
								totalQueryLength);
						closingHits.add(h);
					}

					if (i < frameIDs.length()) {

						// computing offset at query start
						int offset = frameID - lastFrameID;

						// resetting parameters for recording next frame hit
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

			frame = frame_Direction == Frame_Direction.Positiv ? (queryStart % 3) + 1 : -(((queryStart - 1) % 3) + 1);
			queryStart = frame_Direction == Frame_Direction.Positiv ? queryStart - 1 : totalQueryLength - (queryStart);

			// deriving alignment strings
			String[] ali = { aliResult[0], aliResult[1] };
			String[] aliStrings = new CompressAlignment().run(ali);

			// computing edit operations
			Vector<Byte> editOperations = new DAACompressAlignment().run(ali);

			// adding metaInfo
			Object[] metaInfo = { new Integer(queryStart), new Integer(frame), editOperations };

			// assessing alignment properties
			int[] aliStats = new AlignmentEvaluater().run(ali, scoringMatrix);
			int queryLength = aliStats[aliStats.length - 1];

			// generating hit
			Hit h = new Hit(-1, refStart, refEnd, bitScore, rawScore, hit.getFile_pointer(), hit.getAccessPoint(), queryStart, refLength, queryLength,
					hit.getSubjectID());
			h.setFrame(frame);
			h.setHitType(HitType.Synthetic);
			h.setMetaInfo(metaInfo);
			h.copyAliStrings(aliStrings);
			h.setAlignmentStats(aliStats);

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
							if (c == '>') {
								doBreak = true;
								break;
							}
							if (c != '\n')
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

	}

	public class HitRun_Comparator implements Comparator<Hit_Run> {

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

	private String reverseComplementString(String s) {
		String rev = new StringBuffer(s).reverse().toString();
		String revComp = rev.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
		return revComp;
	}

	public void checkRun(Hit_Run run, RandomAccessFile rafSAM, RandomAccessFile rafDAA, String query) {

		String q = run.getFrameDirection() == Frame_Direction.Negativ ? reverseComplementString(query) : query;

		// for (int i = 0; i < 3; i++)
		// System.out.println(new CodonTranslator().translate(q.substring(i)));

		for (Hit h : run.getHitRun()) {
			String[] alignments = h.getAccessPoint() != null ? h.getAlignmentStrings(rafDAA, daaReader) : h.getAlignmentStrings(rafSAM);
			Object[] ali = new ReconstructAlignment(scoringMatrix).run(alignments[1], alignments[0], alignments[2]);
			String qAli = ((String) ali[4]).replaceAll("-", "");
			int start = run.getFrameDirection() == Frame_Direction.Negativ ? q.length() - h.getQuery_start() - 1 : h.getQuery_start();
			String qDNA = q.substring(start, start + h.getQuery_length() * 3);
			String qTrans = new CodonTranslator().translate(qDNA);
			if (!qAli.equals(qTrans)) {
				System.out.println("ERROR:");
				System.out.println(run.getReadID() + " " + run.getGi());
				System.out.println(h);
				System.out.println(qAli + "\n" + qTrans);
				System.out.println("---");
				// printRun(run, rafSAM, rafDAA);
				// System.exit(0);
			}

		}

	}

	public void printRun(Hit_Run run, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		System.out.println("Start: " + run.getHitRun().firstElement().getQuery_start() + " End: " + run.getHitRun().lastElement().getRef_end());
		for (Hit h : run.getHitRun()) {
			h.print("");
			String[] alignments = h.getAccessPoint() != null ? h.getAlignmentStrings(rafDAA, daaReader) : h.getAlignmentStrings(rafSAM);
			Object[] ali = new ReconstructAlignment(scoringMatrix).run(alignments[1], alignments[0], alignments[2]);
			System.out.println(ali[4] + "\n" + ali[5]);
		}
	}

}
