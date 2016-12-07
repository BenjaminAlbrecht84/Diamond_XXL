package pipeline.post.mode_two;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import io.Dmnd_IndexReader;
import io.Fastq_Reader;
import io.SAM_Writer;
import io.daa.DAA_Reader;
import pipeline.post.Alignment_Merger;
import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.HitRun_Writer;
import pipeline.post.HitToSamConverter;
import pipeline.post.Hit_Run;
import pipeline.post.ReadHits;
import pipeline.post.Hit.HitType;
import pipeline.post.mode_one.Alignment_Generator_inParallel;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.SparseString;
import util.ScoringMatrix;
import util.his_algorithm.algorithm.Algorithm_JacobsenVo;

public class Alignment_Generator {

	private CountDownLatch latch;
	private AtomicInteger processedReads = new AtomicInteger(0);
	private int numOfReads, last_p = 0;
	private ConcurrentHashMap<String, Long> readToPointer;
	private Dmnd_IndexReader dmndReader;

	private ConcurrentHashMap<String, ReadHits> readMap;
	private HitRun_Rater runRater;
	private File queryFile, refFile, sam_file;
	private DAA_Reader daaReader;
	private double lambda, k, maxEValue;
	private int minSumScore;
	private int step, length;
	private HitRun_Writer hitRunWriter;
	private HitToSamConverter samConverter;
	private ScoringMatrix matrix;

	public Alignment_Generator(ConcurrentHashMap<String, ReadHits> readMap, HitRun_Rater scorer, File queryFile, File refFile, File sam_file,
			DAA_Reader daaReader, ScoringMatrix matrix, double lambda, double k, int step, int length, HitRun_Writer hitRunWriter,
			HitToSamConverter samConverter, double eValue, int minSumScore) {
		this.readMap = readMap;
		this.runRater = scorer;
		this.queryFile = queryFile;
		this.refFile = refFile;
		this.sam_file = sam_file;
		this.daaReader = daaReader;
		this.matrix = matrix;
		this.lambda = lambda;
		this.k = k;
		this.step = step;
		this.length = length;
		this.hitRunWriter = hitRunWriter;
		this.samConverter = samConverter;
		this.maxEValue = eValue;
		this.minSumScore = minSumScore;
	}

	public void run(int cores) {

		// pre-computing file pointers
		readToPointer = new Fastq_Reader().parseReadIDs(queryFile);
		dmndReader = new Dmnd_IndexReader(refFile);
		dmndReader.createIndex();

		numOfReads = readMap.keySet().size();
		int chunk = (int) Math.ceil((double) numOfReads / (double) cores);

		System.out.println("Processing " + numOfReads + " read-hits...");

		// generating multiple scorer
		Vector<Scorer> allScorer = new Vector<Scorer>();
		Vector<String> chunkfOfReads = new Vector<String>();
		for (String read_name : readMap.keySet()) {
			chunkfOfReads.add(read_name);
			if (!chunkfOfReads.isEmpty() && chunkfOfReads.size() % chunk == 0) {
				allScorer.add(new Scorer(chunkfOfReads));
				chunkfOfReads = new Vector<String>();
			}
		}
		if (!chunkfOfReads.isEmpty()) {
			allScorer.add(new Scorer(chunkfOfReads));
			chunkfOfReads = new Vector<String>();
		}

		// running scorer in parallel
		latch = new CountDownLatch(allScorer.size());
		ExecutorService executor = Executors.newFixedThreadPool(cores);
		for (Scorer thread : allScorer)
			executor.execute(thread);

		// awaiting termination
		try {
			latch.await();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		executor.shutdown();

		System.out.println("OUTPUT>" + 100 + "% (" + processedReads + "/" + numOfReads + ") of the read-hits processed.");

	}

	public class Scorer implements Runnable {

		private final int bufferSize = 10000;
		private Vector<Hit_Run> hitRunBuffer;

		private Vector<String> readNames;
		private Alignment_Completer_Single aliCompleter;

		public Scorer(Vector<String> readNames) {
			this.readNames = readNames;
			this.aliCompleter = new Alignment_Completer_Single(queryFile, refFile, sam_file, daaReader, matrix, lambda, k, runRater, step,
					readToPointer, dmndReader, hitRunWriter, samConverter, maxEValue, minSumScore);
		}

		public void run() {

			hitRunBuffer = new Vector<Hit_Run>();

			try {

				RandomAccessFile rafSAM = new RandomAccessFile(sam_file, "r");
				RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;

				try {

					for (String read_name : readNames) {

						ConcurrentHashMap<SparseString, ConcurrentHashMap<Integer, Vector<Hit>>> hitMap = readMap.get(read_name).getHitMap();
						Vector<SparseString> consideredGIs = new Vector<SparseString>();
						for (SparseString gi : hitMap.keySet()) {

							consideredGIs.add(gi);

							int bestScore = Integer.MIN_VALUE;
							int bestLength = 0;
							int bestRawScore = 0;
							int bestRefLength = 0;
							int bestSumScore = 0;
							double bestEValue = 0;
							Frame_Direction frameDir = Frame_Direction.Positiv;
							Vector<Hit> bestHIS = null;

							// calculate HIS for +&- strand
							for (Frame_Direction dir : Frame_Direction.values()) {

								// collect all hits for the respecting frame
								// direction
								Vector<Hit> allHits = new Vector<Hit>();
								for (int frame : hitMap.get(gi).keySet()) {
									if ((dir == Frame_Direction.Positiv && frame > 0) || (dir == Frame_Direction.Negativ && frame < 0)) {
										Vector<Hit> resFrameHits = mergeOverlaps(hitMap.get(gi).get(frame), dir, runRater, sam_file, matrix, rafSAM,
												rafDAA);
										allHits.addAll(resFrameHits);
									}
								}

								if (!allHits.isEmpty()) {

									// adapting readIDs
									Collections.sort(allHits, new Hit_Comparator());
									HashMap<Hit, Integer> hitToID = new HashMap<Hit, Integer>();
									for (Hit h : allHits)
										hitToID.put(h, h.getId());

									// calculating heaviest increasing subsequence
									Hit[] seq = generateHitSequence(allHits);
									Vector<Hit> his = new Algorithm_JacobsenVo().run(seq, runRater, dir, rafSAM, rafDAA, gi, read_name);
									for (Hit h : his)
										h.setId(hitToID.get(h));

									// store best HIS
									Object[] res = runRater.run(his, dir, rafSAM, rafDAA, read_name, true, gi, read_name);
									if (bestHIS == null || (int) res[2] > bestScore) {
										bestHIS = his;
										bestScore = (int) res[2];
										bestSumScore = (int) res[0];
										bestLength = (int) res[1];
										bestRawScore = (int) res[2];
										bestRefLength = (int) res[3];
										bestEValue = (double) res[4];
										frameDir = dir;
									}

								}

							}

							Vector<Hit> bestHISClone = new Vector<Hit>();
							for (Hit h : bestHIS) {
								Hit hCloned = new Hit(h);
								if (h.getHitType() == HitType.Merged) {
									hCloned.copyAliStrings(h.getAliStrings());
									Object[] metaInfo = { read_name, gi };
									hCloned.setMetaInfo(metaInfo);
								}
								bestHISClone.add(hCloned);
							}

							// adding hit to buffer
							if (bestEValue > maxEValue && bestSumScore > minSumScore) {
								Hit_Run hitRun = new Hit_Run(bestHISClone, new String(read_name), new SparseString(gi), new Integer(bestSumScore),
										new Integer(bestLength), new Integer(bestRawScore), frameDir, new Integer(bestRefLength),
										new Double(bestEValue));
								addHit(hitRun);
							}

							// freeing memory
							for (int frame : hitMap.get(gi).keySet())
								readMap.get(read_name).freeFrameHits(gi, frame);

						}

						// freeing Memory
						for (SparseString gi : consideredGIs)
							readMap.get(read_name).freeGiHits(gi);
						readMap.remove(read_name);

						int p = (int) Math.round(((double) processedReads.incrementAndGet() / (double) numOfReads) * 100.);
						reportProgress(p);

					}

				} finally {
					rafSAM.close();
				}

			} catch (Exception e) {
				e.printStackTrace();
			}

			completeAlignments();
			latch.countDown();

		}

		private Hit[] generateHitSequence(Vector<Hit> allHits) {

			Integer mulKey = null;
			HashMap<Integer, Vector<Hit>> idToHits = new HashMap<Integer, Vector<Hit>>();
			for (Hit h : allHits) {
				if (!idToHits.containsKey(h.getId()))
					idToHits.put(h.getId(), new Vector<Hit>());
				idToHits.get(h.getId()).add(h);
				if (idToHits.get(h.getId()).size() > 1)
					mulKey = h.getId();
			}

			if (mulKey == null)
				return allHits.toArray(new Hit[allHits.size()]);

			Vector<Hit> mulHits = idToHits.get(mulKey);
			int counter = mulHits.firstElement().getId();
			for (Hit h : allHits) {
				if (mulHits.contains(h))
					h.setId(counter++);
				else if (h.getId() > mulKey)
					h.setId(h.getId() + mulHits.size());
			}

			return generateHitSequence(allHits);

		}

		private void addHit(Hit_Run hitRun) {
			hitRunBuffer.add(hitRun);
			if (hitRunBuffer.size() % bufferSize == 0)
				completeAlignments();
		}

		private void completeAlignments() {
			aliCompleter.run(hitRunBuffer);

			// freeing Memory
			for (Hit_Run run : hitRunBuffer) {
				for (Hit h : run.getHitRun())
					h.freeMemory();
			}

			// resetting buffer
			hitRunBuffer = new Vector<Hit_Run>();
		}

	}

	private synchronized void reportProgress(int p) {
		if (p != 100 && p != last_p && p % 10 == 0) {
			System.out.println("OUTPUT>" + p + "% (" + processedReads + "/" + numOfReads + ") of the read-hits processed.");
			last_p = p;
		}
	}

	private Vector<Hit> mergeOverlaps(Vector<Hit> frameHits, Frame_Direction dir, HitRun_Rater scorer, File sam_file, ScoringMatrix matrix,
			RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		Collections.sort(frameHits, new Hit_Comparator());

		// System.out.println("\n" + dir);
		// for (Hit h : frameHits) {
		// try {
		// RandomAccessFile raf = new RandomAccessFile(sam_file, "r");
		// raf.seek(h.getFile_pointer());
		// String line = raf.readLine();
		// System.out.println(h.getRef_start());
		// System.out.println(line);
		// raf.close();
		// } catch (Exception e) {
		//
		// }
		// }

		Vector<Hit> resHits = new Vector<Hit>();
		resHits.add(frameHits.get(0));
		for (int i = 1; i < frameHits.size(); i++) {
			Hit hL = resHits.lastElement();
			Hit hR = frameHits.get(i);
			if (hL.getRef_end() >= hR.getRef_start()) {

				if (hL.getRef_start() < hR.getRef_start() || hL.getQuery_length() < hR.getQuery_length()) {

					Hit hM = new Alignment_Merger(scorer, sam_file, daaReader, step, length).mergeTwoHits(hL, hR, matrix, rafSAM, rafDAA, null, null,
							null, dir);

					resHits.remove(hL);
					resHits.add(hM);

				}

			} else
				resHits.add(hR);
		}

		return resHits;
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

}
