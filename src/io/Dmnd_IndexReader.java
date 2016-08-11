package io;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

public class Dmnd_IndexReader {

	private File dmndFile;

	// private TreeMap<Integer, Long> giIndex;
	private ConcurrentHashMap<Integer, Long> giIndex;
	private Vector<long[]> seqLocations;

	public Dmnd_IndexReader(File dmndFile) {
		this.dmndFile = dmndFile;
	}

	public void createIndex() {

		parseSeqLocations();
		giIndex = new ConcurrentHashMap<Integer, Long>();
		giIndex.putAll(mapGIs());
		seqLocations = null;

	}

	public HashMap<Integer, Long> getGILocations(Vector<Integer> gIs) {

		HashMap<Integer, Long> gILocations = new HashMap<Integer, Long>();
		for (int gi : gIs)
			gILocations.put(gi, giIndex.get(gi));

		return gILocations;
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

				while ((readChars = is.read(buffer.array())) != -1) {

					for (int i = 0; i < readChars - 15; i += 16) {

						long start = (long) buffer.getLong(i);
						long length = (long) buffer.getLong(i + 8);
						if (length != 0) {
							long[] loc = { start, length };
							seqLocations.add(loc);
						}

					}

				}

			} finally {
				is.close();
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private ConcurrentHashMap<Integer, Long> mapGIs() {

		System.out.println("Indexing reference database...");

		ConcurrentHashMap<Integer, Long> giToPointer = new ConcurrentHashMap<Integer, Long>();
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
							if (val == 124) {
								int gi = Integer.parseInt(giBuffer.toString());
								giBuffer = new StringBuffer();
								giToPointer.put(gi, seqLocations.get(lastK)[0]);
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
						long giStart = (loc[0] + loc[1] + 5) - (allocSize * iterations);

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
							if (val == 124) {
								try {
									int gi = Integer.parseInt(giBuffer.toString());
									giBuffer = new StringBuffer();
									giToPointer.put(gi, loc[0]);
									break;
								} catch (Exception e) {
									e.printStackTrace();
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

		int proc = (int) Math.ceil(((double) (giToPointer.keySet().size()) / (double) seqLocations.size()) * 100.);
		System.out.println("OUTPUT>Finally, " + proc + "% (" + giToPointer.keySet().size() + "/" + seqLocations.size() + ") of the GI's mapped...");

		return giToPointer;

	}

}
