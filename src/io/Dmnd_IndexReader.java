package io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import util.SparseString;

public class Dmnd_IndexReader {

	private File dmndFile;

	private ConcurrentHashMap<SparseString, SparseLong> giIndex;
	private Vector<long[]> seqLocations;

	private boolean giveWarnings = false;
	private int indexChunk = 0, totalIndexChunks = 1, totalNumberOfSeqIDs = 0;

	public Dmnd_IndexReader(File dmndFile) {
		this.dmndFile = dmndFile;
	}

	public Dmnd_IndexReader(File dmndFile, int indexChunk, int totalIndexChunks) {
		this.dmndFile = dmndFile;
		this.indexChunk = indexChunk;
		this.totalIndexChunks = totalIndexChunks;
	}

	public void createIndex() {

		System.out.println("STEP_5>Indexing reference database...");
		long time = System.currentTimeMillis();

		parseSeqLocations();
		giIndex = new ConcurrentHashMap<SparseString, SparseLong>();
		giIndex.putAll(mapGIs());

		long runtime = (System.currentTimeMillis() - time) / 1000;
		int proc = (int) Math.ceil(((double) (giIndex.keySet().size()) / (double) seqLocations.size()) * 100.);
		System.out.println(
				"OUTPUT>" + giIndex.keySet().size() + "/" + totalNumberOfSeqIDs + " SeqIDs mapped. [" + runtime + "s]\n");

		seqLocations = null;

	}

	public HashMap<SparseString, Long> getGILocations(Vector<SparseString> gIs) {
		HashMap<SparseString, Long> gILocations = new HashMap<SparseString, Long>();
		for (SparseString gi : gIs) {
			if (giIndex.containsKey(gi))
				gILocations.put(gi, new Long(giIndex.get(gi).getValue()));
		}
		return gILocations;
	}

	public Long getGILocation(SparseString gi) {
		if (giIndex.containsKey(gi))
			return giIndex.get(gi).getValue();
		return null;
	}

	private void parseSeqLocations() {

		seqLocations = new Vector<long[]>();

		try {

			// parsing sequence bounds
			InputStream is = new BufferedInputStream(new FileInputStream(dmndFile));
			long offset;
			try {
				ByteBuffer buffer = ByteBuffer.allocate(40);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				is.read(buffer.array());
				int id = buffer.getInt(0);
				int build = buffer.getShort(8);
				int vers = buffer.getShort(12);
				int seq = buffer.getInt(16);
				int letters = buffer.getInt(24);
				offset = (long) buffer.getLong(32);

			} finally {
				is.close();
			}

			is = new BufferedInputStream(new FileInputStream(dmndFile));
			try {

				long totalBytes = dmndFile.length() - offset;
				int allocSize = 1024 * 1024;
				ByteBuffer buffer = ByteBuffer.allocate(allocSize);
				buffer.order(ByteOrder.LITTLE_ENDIAN);

				long skipped = is.skip(offset);
				int readChars = 0, lastProc = 0;
				long totalReadBytes = 0;

				if (skipped != offset)
					throw new IllegalArgumentException("ERROR: Too less bytes have been skipped! " + skipped + " " + offset);

				Vector<long[]> allSeqLocations = new Vector<long[]>();
				while ((readChars = is.read(buffer.array())) != -1) {
					for (int i = 0; i < readChars - 15; i += 16) {
						long start = (long) buffer.getLong(i);
						long length = (long) buffer.getLong(i + 8);
						if (length != 0) {
							long[] loc = { start, length };
							allSeqLocations.add(loc);
						}

					}

				}
				totalNumberOfSeqIDs = allSeqLocations.size();

				// storing relevant seqLocations
				int chunkSize = (int) Math.ceil((double) allSeqLocations.size() / (double) totalIndexChunks);
				int start = indexChunk * chunkSize, end = start + chunkSize;
				for (int i = start; i < end; i++) {
					if (i == allSeqLocations.size())
						break;
					seqLocations.addElement(allSeqLocations.get(i));
				}

			} finally {
				is.close();
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private ConcurrentHashMap<SparseString, SparseLong> mapGIs() {

		ConcurrentHashMap<SparseString, SparseLong> giToPointer = new ConcurrentHashMap<SparseString, SparseLong>();
		try {
			InputStream is = new BufferedInputStream(new FileInputStream(dmndFile));
			try {

				int lastProc = 0, checkBound = 100000;
				long totalReadBytes = 0;

				int allocSize = 1024 * 1024;
				long iterations = -1;
				ByteBuffer buffer = ByteBuffer.allocate(allocSize);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				int readBytes = 0, lastK = 0;
				long totalBytes = dmndFile.length();
				StringBuffer giBuffer = new StringBuffer("");
				while ((readBytes = is.read(buffer.array())) != -1) {

					iterations++;

					// complete overlapping giNumber from last buffer
					if (giBuffer.length() != 0) {
						for (int i = 0; i < readBytes; i++) {
							int val = (int) buffer.get(i);
							if (val == 32 || val == 0) {
								SparseString gi = new SparseString(giBuffer.toString());
								giToPointer.put(gi, new SparseLong(seqLocations.get(lastK)[0]));
								giBuffer = new StringBuffer();
								lastK++;
								break;
							}
							giBuffer = giBuffer.append((char) val);
						}
					}

					// find new giNumbers in buffer
					for (int k = lastK; k < seqLocations.size(); k++) {

						// k = seqLocations.get(k) == null ? k + 1 : k;
						long[] loc = seqLocations.get(k);
						long giStart = (loc[0] + loc[1] + 2) - (allocSize * iterations);

						// checking if GI is in current buffer
						if (giStart >= readBytes) {
							lastK = k;
							break;
						}

						boolean doBreak = false;
						for (long i = giStart; i <= readBytes; i++) {

							if (doBreak = (i == readBytes))
								break;

							int val = (int) buffer.get((int) i);
							if (val == 32 || val == 0) {
								try {
									SparseString gi = new SparseString(giBuffer.toString());
									giToPointer.put(gi, new SparseLong(loc[0]));
									giBuffer = new StringBuffer();
									break;
								} catch (Exception e) {
									e.printStackTrace();
									System.exit(0);
								}
							}
							giBuffer = giBuffer.append((char) val);
						}

						lastK = k;
						if (doBreak)
							break;

					}

					if (giBuffer.length() == 0 && lastK == seqLocations.size() - 1)
						break;

					totalReadBytes += readBytes;
					int proc = (int) Math.floor(((double) (totalReadBytes) / (double) totalBytes) * 100.);
					if (proc % 10 == 0 && proc > lastProc) {
						lastProc = proc;
						System.out.println("OUTPUT>" + proc + "% of the index constructed...");
					}

				}

			} finally {
				is.close();
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return giToPointer;

	}

	public class SparseLong {

		private final long l;

		public SparseLong(long l) {
			this.l = l;
		}

		public long getValue() {
			return l;
		}

	}

	public int getNumberOfSequences() {
		try {
			InputStream is = new BufferedInputStream(new FileInputStream(dmndFile));
			try {
				ByteBuffer buffer = ByteBuffer.allocate(40);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				is.read(buffer.array());
				return buffer.getInt(16);
			} finally {
				is.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return -1;
	}

	public ConcurrentHashMap<SparseString, SparseLong> getGiIndex() {
		return giIndex;
	}

}
