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
import io.debug.RejectedWriter;
import pipeline.post.Alignment_Merger;
import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.Hit_Run;
import pipeline.post.ReadHits;
import pipeline.post.Hit.HitType;
import util.HitLine;
import util.ReconstructAlignment;
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
	private RejectedWriter rejectedWriter;

	private ConcurrentHashMap<String, ReadHits> readMap;
	private HitRun_Rater scorer;
	private File sam_file, rej_file_1;
	private DAA_Reader daaReader;
	private ScoringMatrix matrix;
	private double maxEValue;
	private int minSumScore;
	private int step, length;
	private boolean useFilters;

	public Alignment_Generator_inParallel(ConcurrentHashMap<String, ReadHits> readMap, HitRun_Rater scorer, File sam_file, DAA_Reader daaReader,
			ScoringMatrix matrix, double maxEValue, int minSumScore, boolean useFilters, int step, int length, File rej_file_1) {
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
		this.rej_file_1 = rej_file_1;
	}

	public Vector<Hit_Run> run(int cores) {

		Runtime.getRuntime().availableProcessors();

		rejectedWriter = rej_file_1 != null ? new RejectedWriter(rej_file_1) : null;
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
			Vector<Hit_Run> rejectedHitRuns = new Vector<Hit_Run>();
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
									Vector<Hit> his = new Algorithm_JacobsenVo().run(seq, scorer, dir, rafSAM, rafDAA, gi, read_name);
									for (Hit h : his)
										h.setId(hitToID.get(h));

									// filtering his after distance to diagonal line
									if (his.size() > 1 && useFilters) {

										// again computing heaviest increasing subsequence on FILTERED hits
										// Vector<Hit> filteredHits = filterHits(allHits, dir, read_name, gi);
										Vector<Hit> filteredHits = filterHits_AllPairs(allHits, dir, rafSAM, rafDAA, read_name, gi, read_name);

										seq = generateHitSequence(filteredHits);
										his = new Algorithm_JacobsenVo().run(seq, scorer, dir, rafSAM, rafDAA, gi, read_name);
										for (Hit h : his)
											h.setId(hitToID.get(h));

									}

									// storing best result
									Object[] res = scorer.run(his, dir, rafSAM, rafDAA, read_name, true, gi, read_name);
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
										hCloned.copyAliStrings(h.getAliStrings());
										Object[] metaInfo = { read_name, gi };
										hCloned.setMetaInfo(metaInfo);
										mergedHits.add(hCloned);
									}
									bestHISClone.add(hCloned);
								}

								if (!useFilters || bestSumScore > minSumScore) {
									Hit_Run hitRun = new Hit_Run(bestHISClone, new String(read_name), new SparseString(gi), new Integer(bestSumScore),
											new Integer(bestLength), new Integer(bestRawScore), frameDir, new Integer(bestRefLength),
											new Double(bestEValue));
									localHitRuns.add(hitRun);
								} else if (rejectedWriter != null) {
									Hit_Run hitRun = new Hit_Run(bestHISClone, new String(read_name), new SparseString(gi), new Integer(bestSumScore),
											new Integer(bestLength), new Integer(bestRawScore), frameDir, new Integer(bestRefLength),
											new Double(bestEValue));
									rejectedHitRuns.add(hitRun);
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

			// reporting rejected runs
			if (rejectedWriter != null)
				rejectedWriter.run(rejectedHitRuns);

			// writing new hits to file
			samWriter.run(mergedHits);

			latch.countDown();

		}

		private Vector<Hit> filterHits_AllPairs(Vector<Hit> hits, Frame_Direction dir, RandomAccessFile rafSAM, RandomAccessFile rafDAA,
				String readID, SparseString gi, String read_name) {

			if (hits.size() < 2)
				return hits;

			HitLine bestHitLine = null;
			int bestScore = Integer.MIN_VALUE;
			Vector<Hit> bestLinearHits = new Vector<Hit>();
			for (int i = 0; i < hits.size() - 1; i++) {
				for (int j = i + 1; j < hits.size(); j++) {

					Hit h1 = hits.get(i);
					Hit h2 = hits.get(j);

					// fitting 2D line with fixed slope of 45 degrees
					HitLine hitLine = new HitLine(step, h1, h2, dir);

					// filtering method follows the alignment refinement step of GraphMap
					// -> see GraphMap-Paper: "Fast and sensitive mapping of nanopore sequencing reads with GraphMap"

					double maxDist = 0.45 * (double) hits.firstElement().getRef_length() * Math.sqrt(2.) / 2.;

					double sumDist = 0, n = 0;
					for (Hit h : hits) {
						if (hitLine.getDistance(h, dir) <= maxDist) {
							n++;
							sumDist += hitLine.getDistance(h, dir);
						}
					}

					double c = 3. * (sumDist / n);
					Vector<Hit> filteredHits = new Vector<Hit>();
					for (Hit h : hits) {
						if (hitLine.getDistance(h, dir) <= c)
							filteredHits.add(h);
					}

					if (filteredHits.isEmpty()) {
						if (h1.getBitScore() > h2.getBitScore())
							filteredHits.add(h1);
						else
							filteredHits.add(h2);
					}

					Object[] res = scorer.run(filteredHits, dir, rafSAM, rafDAA, read_name, true, gi, read_name);
					int score = (int) res[0];
					if (score > bestScore) {
						bestHitLine = hitLine;
						bestScore = score;
						bestLinearHits = filteredHits;
					}

				}
			}

			// if (bestLinearHits.size() < hits.size()) {
			// System.out.println(read_name + " " + gi);
			// System.out.println(bestHitLine.getLog());
			// System.out.println("SIZE: " + bestLinearHits.size());
			// for (Hit h : hits)
			// System.out.println(bestLinearHits.contains(h) + " " + h);
			// }

			return bestLinearHits;

		}

		private Vector<Hit> filterHits(Vector<Hit> hits, Frame_Direction dir, String read_name, SparseString gi) {

			if (hits.size() < 2)
				return hits;

			// fitting 2D line with fixed slope of 45 degrees
			HitLine hitLine = new HitLine(step, hits, dir);

			// filtering method follows the alignment refinement step of GraphMap
			// -> see GraphMap-Paper: "Fast and sensitive mapping of nanopore sequencing reads with GraphMap"

			double maxDist = 0.45 * (double) hits.firstElement().getRef_length() * Math.sqrt(2.) / 2.;

			double sumDist = 0, n = 0;
			for (Hit h : hits) {
				if (hitLine.getDistance(h, dir) <= maxDist) {
					n++;
					sumDist += hitLine.getDistance(h, dir);
				}
			}

			double c = 9. * (sumDist / n);
			Vector<Hit> filteredHits = new Vector<Hit>();
			for (Hit h : hits) {
				if (hitLine.getDistance(h, dir) <= c)
					filteredHits.add(h);
			}

			// if (filteredHits.size() < hits.size()) {
			// System.out.println(read_name + " " + gi);
			// System.out.println(hitLine.getLog());
			// System.out.println("SIZE: " + filteredHits.size());
			// for (Hit h : hits)
			// System.out.println(filteredHits.contains(h) + " " + h);
			// }

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

				if (hL.getRef_start() <= hR.getRef_start() || hL.getQuery_length() < hR.getQuery_length()) {

					Hit hM = new Alignment_Merger(scorer, sam_file, daaReader, step, length).mergeTwoHits(hL, hR, matrix, rafSAM, rafDAA, null, null,
							null, dir);

					int offset_L = Math.abs(hL.getId() * step);
					int qStart_L = dir == Frame_Direction.Positiv ? hL.getQuery_start() + offset_L - 1 : length - hL.getQuery_start() - offset_L;
					int qEnd_L = qStart_L + hL.getQuery_length() * 3;

					int offset_R = Math.abs(hR.getId() * step);
					int qStart_R = dir == Frame_Direction.Positiv ? hR.getQuery_start() + offset_R - 1 : length - hR.getQuery_start() - offset_R;
					int qEnd_R = qStart_R + hR.getQuery_length() * 3;

					int offset_M = Math.abs(hM.getId() * step);
					int qStart_M = dir == Frame_Direction.Positiv ? hM.getQuery_start() + offset_M - 1 : length - hM.getQuery_start() - offset_M;
					int qEnd_M = qStart_M + hM.getQuery_length() * 3;

					int qEnd = qEnd_R > qEnd_L ? qEnd_R : qEnd_L;
					if (qEnd_M == qEnd) {

						resHits.remove(hL);
						resHits.add(hM);

						// System.out.println();
						// printHit(hL, rafSAM, rafDAA);
						// printHit(hR, rafSAM, rafDAA);
						// printHit(hM, rafSAM, rafDAA);

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

	public void printHit(Hit h, RandomAccessFile rafSAM, RandomAccessFile rafDAA) {
		h.print("");

		String[] alignments = h.getAccessPoint() != null ? h.getAlignmentStrings(rafDAA, daaReader) : h.getAlignmentStrings(rafSAM);
		Object[] ali = new ReconstructAlignment(matrix).run(alignments[1], alignments[0], alignments[2]);
		System.out.println(ali[4] + "\n" + ali[5]);
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
