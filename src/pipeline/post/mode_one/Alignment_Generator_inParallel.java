package pipeline.post.mode_one;

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

import io.SAM_Writer;
import io.daa.DAA_Reader;
import pipeline.post.Alignment_Merger;
import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.Hit_Run;
import pipeline.post.ReadHits;
import pipeline.post.Hit.HitType;
import util.HitLine;
import util.SparseString;
import util.ScoringMatrix;
import util.his_algorithm.algorithm.Algorithm_JacobsenVo;

public class Alignment_Generator_inParallel {

	public enum Frame_Direction {
		Positiv, Negativ
	}

	private Vector<Hit_Run> hitRuns;
	private CountDownLatch latch;
	private AtomicInteger processedReads = new AtomicInteger(0);
	private int numOfReads, last_p = 0;
	private SAM_Writer samWriter;

	private ConcurrentHashMap<String, ReadHits> readMap;
	private HitRun_Rater scorer;
	private File sam_file;
	private DAA_Reader daaReader;
	private ScoringMatrix matrix;
	private double maxEValue;
	private int minSumScore;
	private int step, length;
	private boolean useFilters;

	public Alignment_Generator_inParallel(ConcurrentHashMap<String, ReadHits> readMap, HitRun_Rater scorer, File sam_file, DAA_Reader daaReader,
			ScoringMatrix matrix, double maxEValue, int minSumScore, boolean useFilters, int step, int length) {
		this.readMap = readMap;
		this.scorer = scorer;
		this.sam_file = sam_file;
		this.daaReader = daaReader;
		this.matrix = matrix;
		this.maxEValue = maxEValue;
		this.minSumScore = minSumScore;
		this.useFilters = useFilters;
		this.step = step;
		this.length = length;
	}

	public Vector<Hit_Run> run(int cores) {

		Runtime.getRuntime().availableProcessors();

		samWriter = new SAM_Writer(sam_file, daaReader);
		hitRuns = new Vector<Hit_Run>();
		numOfReads = readMap.keySet().size();
		int chunk = (int) Math.ceil((double) numOfReads / (double) cores);
		chunk = chunk > 100 ? 100 : chunk;

		System.out.println("STEP_4>Processing " + numOfReads + " read-hit(s)...");
		long time = System.currentTimeMillis();

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

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.out.println("OUTPUT>" + 100 + "% (" + processedReads + "/" + numOfReads + ") of the read-hits processed. [" + runtime + "s]\n");

		return hitRuns;

	}

	private synchronized void reportLocalHitRuns(Hit_Run run) {
		hitRuns.add(run);
	}

	public class Scorer implements Runnable {

		private Vector<String> readNames;

		public Scorer(Vector<String> readNames) {
			this.readNames = readNames;
		}

		public void run() {

			Vector<Hit_Run> localHitRuns = new Vector<Hit_Run>();
			ConcurrentLinkedQueue<Hit> mergedHits = new ConcurrentLinkedQueue<Hit>();

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

							// calculating HIS for POSITIV and NEGATIVE strand
							for (Frame_Direction dir : Frame_Direction.values()) {

								// collecting all hits for the respecting frame direction
								Vector<Hit> allHits = new Vector<Hit>();
								for (int frame : hitMap.get(gi).keySet()) {
									if ((dir == Frame_Direction.Positiv && frame > 0) || (dir == Frame_Direction.Negativ && frame < 0)) {
										Vector<Hit> resFrameHits = mergeOverlaps(hitMap.get(gi).get(frame), dir, scorer, sam_file, matrix, rafSAM,
												rafDAA);
										allHits.addAll(resFrameHits);
									}
								}

								if (!allHits.isEmpty()) {

									Collections.sort(allHits, new Hit_Comparator());
									HashMap<Hit, Integer> hitToID = new HashMap<Hit, Integer>();
									for (Hit h : allHits)
										hitToID.put(h, h.getId());

									// computing heaviest increasing subsequence on ALL hits
									Hit[] seq = generateHitSequence(allHits);
									Vector<Hit> his = new Algorithm_JacobsenVo().run(seq, scorer, dir, rafSAM, rafDAA, read_name);
									for (Hit h : his)
										h.setId(hitToID.get(h));

									// filtering his after distance to diagonal line
									if (useFilters) {

										// fitting 2D line with fixed slope of 45 degrees
										HitLine hitLine = new HitLine(step, his);

										// again computing heaviest increasing subsequence on FILTERED hits
										Vector<Hit> filteredHits = filterHits(hitLine, allHits);
										seq = generateHitSequence(filteredHits);
										his = new Algorithm_JacobsenVo().run(seq, scorer, dir, rafSAM, rafDAA, read_name);
										for (Hit h : his)
											h.setId(hitToID.get(h));

									}

									// storing best result
									Object[] res = scorer.run(his, dir, rafSAM, rafDAA, read_name, true);
									if (!his.isEmpty() && (bestHIS == null || (int) res[2] > bestScore)) {
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

							// storing best HIS
							if (bestHIS != null && !bestHIS.isEmpty()) {

								Vector<Hit> bestHISClone = new Vector<Hit>();
								for (Hit h : bestHIS) {
									Hit hCloned = new Hit(h);
									if (h.getHitType() == HitType.Merged) {
										// hCloned.setHitType(HitType.Merged);
										hCloned.copyAliStrings(h.getAliStrings());
										Object[] metaInfo = { read_name, gi };
										hCloned.setMetaInfo(metaInfo);
										mergedHits.add(hCloned);
									}
									bestHISClone.add(hCloned);
								}

								if (!useFilters || (bestEValue < maxEValue && bestSumScore > minSumScore)) {
									Hit_Run hitRun = new Hit_Run(bestHISClone, new String(read_name), new SparseString(gi), new Integer(bestSumScore),
											new Integer(bestLength), new Integer(bestRawScore), frameDir, new Integer(bestRefLength),
											new Double(bestEValue));
									localHitRuns.add(hitRun);
								}

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
					if (rafDAA != null)
						rafDAA.close();
				}

			} catch (Exception e) {
				e.printStackTrace();
			}

			// reporting best runs
			for (Hit_Run run : localHitRuns)
				reportLocalHitRuns(run);

			// writing new hits to file
			samWriter.run(mergedHits);

			latch.countDown();

		}

		private Vector<Hit> filterHits(HitLine hitLine, Vector<Hit> hits) {

			// filtering method follows the alignment refinement step of GraphMap
			// -> see GraphMap-Paper: "Fast and sensitive mapping of nanopore sequencing reads with GraphMap"

			double maxDist = 0.45 * (double) hits.firstElement().getRef_length() * Math.sqrt(2.) / 2.;

			double sumDist = 0, n = 0;
			for (Hit h : hits) {
				if (hitLine.getDistance(h) <= maxDist) {
					n++;
					sumDist += hitLine.getDistance(h);
				}
			}

			double c = 3. * (sumDist / n);
			Vector<Hit> filteredHits = new Vector<Hit>();
			for (Hit h : hits) {
				if (hitLine.getDistance(h) <= c)
					filteredHits.add(h);
			}

			return filteredHits;
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

		Vector<Hit> recHits = new Vector<Hit>();
		Vector<Hit> resHits = new Vector<Hit>();
		resHits.add(frameHits.get(0));
		for (int i = 1; i < frameHits.size(); i++) {

			Hit hL = resHits.lastElement();
			Hit hR = frameHits.get(i);

			if (hL.getRef_end() >= hR.getRef_start()) {

				if (hL.getRef_start() < hR.getRef_start() || hL.getQuery_length() < hR.getQuery_length()) {

					Hit hM = new Alignment_Merger(scorer, sam_file, daaReader).mergeTwoHits(hL, hR, matrix, rafSAM, rafDAA, null, null, null);

					int offset_R = Math.abs(hR.getId() * step);
					int qStart_R = dir == Frame_Direction.Positiv ? hR.getQuery_start() + offset_R - 1 : length - hR.getQuery_start() - offset_R;
					int qEnd_R = qStart_R + hR.getQuery_length() * 3;

					int offset_M = Math.abs(hM.getId() * step);
					int qStart_M = dir == Frame_Direction.Positiv ? hM.getQuery_start() + offset_M - 1 : length - hM.getQuery_start() - offset_M;
					int qEnd_M = qStart_M + hM.getQuery_length() * 3;

					if (qEnd_M == qEnd_R) {
						resHits.remove(hL);
						resHits.add(hM);
					} else
						recHits.add(hR);

				}

			} else
				resHits.add(hR);
		}

		if (!recHits.isEmpty())
			resHits.addAll(mergeOverlaps(recHits, dir, scorer, sam_file, matrix, rafSAM, rafDAA));

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
