package io.daa;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicLong;

import pipeline.post.Hit;
import pipeline.post.Hit_Run;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;

public class DAA_Writer {

	private AtomicLong queryRecords = new AtomicLong(0), aliBlockSize = new AtomicLong(0);
	private HashMap<Character, Integer> nucToIndex;

	private File out, samFile;
	private DAA_Reader daaReader;

	public DAA_Writer(File out, File samFile, DAA_Reader daaReader) {
		this.out = out;
		this.samFile = samFile;
		this.daaReader = daaReader;

		nucToIndex = new HashMap<Character, Integer>();
		nucToIndex.put('A', 0);
		nucToIndex.put('C', 1);
		nucToIndex.put('G', 2);
		nucToIndex.put('T', 3);

		writeHeader();
	}

	private void writeHeader() {
		try {
			InputStream is = new BufferedInputStream(new FileInputStream(daaReader.getDAAFile()));
			try {

				long headerSize = daaReader.getDAAHeader().getHeaderSize();
				ByteBuffer buffer = ByteBuffer.allocate((int) headerSize);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				is.read(buffer.array());
				writeInFile(buffer.array(), buffer.capacity(), false);

				// System.out.println("Header: " + out.length());

			} finally {
				is.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void run(Vector<Hit_Run> runs, HashMap<String, String> readIDToSeq) {

		try {

			ArrayList<Byte> byteBuffer = new ArrayList<Byte>();

			RandomAccessFile rafSAM = new RandomAccessFile(samFile, "r");
			RandomAccessFile rafDAA = daaReader != null ? new RandomAccessFile(daaReader.getDAAFile(), "r") : null;

			HashMap<String, byte[]> readIDToPackedSeq = new HashMap<String, byte[]>();

			try {

				Vector<int[]> allocPairs = new Vector<int[]>();
				Collections.sort(runs, new HitRunComparator());
				String lastReadID = "";
				int begin = -1, alloc = -1;

				for (int i = 0; i < runs.size(); i++) {

					Hit_Run run = runs.get(i);

					if (lastReadID.isEmpty() || !run.getReadID().equals(lastReadID)) {

						if (alloc != -1) {
							alloc = byteBuffer.size() - 4 - begin;
							int[] allocPair = { begin, alloc };
							allocPairs.add(allocPair);
						}

						begin = byteBuffer.size();

						alloc = 0;
						write(byteBuffer, readLittleEndian(alloc));

						int totalQueryLength = readIDToSeq.get(run.getReadID()).length();
						write(byteBuffer, readLittleEndian(totalQueryLength));

						String queryName = run.getReadID();
						write(byteBuffer, readLittleEndian(queryName));
						write(byteBuffer, (byte) 0);

						byte nFlag = 0;
						write(byteBuffer, nFlag);

						if (!readIDToPackedSeq.containsKey(queryName)) {
							byte[] packedSequence = packSequence(readIDToSeq.get(queryName));
							readIDToPackedSeq.put(queryName, packedSequence);
						}
						byte[] packedSequence = readIDToPackedSeq.get(run.getReadID());
						write(byteBuffer, packedSequence);

					}

					for (Hit h : run.getHitRun()) {

						int subjectID = h.getSubjectID();
						write(byteBuffer, readLittleEndian(subjectID));

						byte typeFlags = 1 << 1;
						typeFlags |= 1 << 3;
						typeFlags |= 1 << 5;
						typeFlags |= run.getFrameDirection() == Frame_Direction.Positiv ? 0 : 1 << 6;
						write(byteBuffer, typeFlags);

						int rawScore = h.getRawScore();
						write(byteBuffer, readLittleEndian(rawScore));

						int queryStart = h.getQuery_start() - 1;
						write(byteBuffer, readLittleEndian(queryStart));

						int refStart = h.getRef_start();
						write(byteBuffer, readLittleEndian(refStart));

						Vector<Byte> editOperations = (Vector<Byte>) h.getMetaInfo()[2];
						for (byte op : editOperations)
							write(byteBuffer, op);
						write(byteBuffer, (byte) 0);

					}

					lastReadID = run.getReadID();

				}

				if (alloc != -1) {

					alloc = byteBuffer.size() - 4 - begin;
					int[] allocPair = { begin, alloc };
					allocPairs.add(allocPair);

					byte[] stream = new byte[byteBuffer.size()];
					for (int i = 0; i < byteBuffer.size(); i++)
						stream[i] = byteBuffer.get(i);

					for (int[] p : allocPairs) {
						int counter = 0;
						for (byte b : readLittleEndian(p[1])) {
							int pos = p[0] + (counter++);
							stream[pos] = b;
						}
					}

					writeInFile(stream, stream.length, true);

					aliBlockSize.getAndAdd(stream.length);
					queryRecords.getAndAdd(allocPairs.size());

				}

			} finally {
				rafSAM.close();
				if (rafDAA != null)
					rafDAA.close();
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void write(ArrayList<Byte> byteBuffer, byte b) {
		byteBuffer.add(b);
	}

	private void write(ArrayList<Byte> byteBuffer, byte[] bytes) {
		for (byte b : bytes)
			byteBuffer.add(b);
	}

	private byte[] readLittleEndian(String s) {
		ByteBuffer buffer = ByteBuffer.allocate(s.length());
		buffer.order(ByteOrder.LITTLE_ENDIAN);
		for (int i = 0; i < s.length(); i++)
			buffer.put((byte) s.charAt(i));
		return buffer.array();
	}

	private byte[] readLittleEndian(int o) {
		ByteBuffer buffer = ByteBuffer.allocate(4);
		buffer.order(ByteOrder.LITTLE_ENDIAN);
		buffer.putInt(o);
		return buffer.array();
	}

	private byte[] readLittleEndian(long o) {
		ByteBuffer buffer = ByteBuffer.allocate(8);
		buffer.order(ByteOrder.LITTLE_ENDIAN);
		buffer.putLong(o);
		return buffer.array();
	}

	private byte[] packSequence(String dna) {

		Vector<Byte> packed = new Vector<Byte>();
		byte p = 0;
		for (int i = 0; i < dna.length(); i++) {
			char c = dna.charAt(i);
			byte b = 0;
			b |= nucToIndex.get(c) << (i * 2) % 8;
			p |= b;
			if (i == dna.length() - 1 || (((i + 1) * 2) % 8 == 0 && i != 0)) {
				packed.add(p);
				p = 0;
			}
		}

		byte[] a = new byte[packed.size()];
		for (int i = 0; i < packed.size(); i++)
			a[i] = packed.get(i);

		return a;

	}

	public void finish() {
		try {

			// finishing alignment block
			writeInFile(readLittleEndian((int) 0), 4, true);
			aliBlockSize.getAndAdd(4);

			// inserting reference names and reference lengths
			InputStream is = new BufferedInputStream(new FileInputStream(daaReader.getDAAFile()));
			ByteBuffer buffer = ByteBuffer.allocate(2448);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			is.read(buffer.array());
			long blockSize = buffer.getLong(144);
			is.close();

			RandomAccessFile raf = new RandomAccessFile(daaReader.getDAAFile(), "r");
			try {
				long headerSize = daaReader.getDAAHeader().getHeaderSize();
				long offset = headerSize + blockSize;
				raf.seek(offset);
				int bufferSize = 1024;
				buffer = ByteBuffer.allocate(bufferSize);
				buffer.order(ByteOrder.LITTLE_ENDIAN);
				int readChars = 0;
				while ((readChars = raf.read(buffer.array())) != -1)
					writeInFile(buffer.array(), readChars, true);
			} finally {
				raf.close();
			}

			// updating query records
			writeByteInFile(readLittleEndian(queryRecords.get()), 56);

			// updating alignment block size
			writeByteInFile(readLittleEndian(aliBlockSize.get()), 144);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private synchronized void writeByteInFile(byte[] b, long pos) {
		try {
			RandomAccessFile raf = new RandomAccessFile(out, "rw");
			try {
				raf.seek(pos);
				raf.write(b);
			} finally {
				raf.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private synchronized void writeInFile(byte[] b, int len, boolean append) {
		try {
			OutputStream output = null;
			try {
				output = new BufferedOutputStream(new FileOutputStream(out, append));
				output.write(b, 0, len);
			} finally {
				output.close();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public class HitRunComparator implements Comparator<Hit_Run> {

		@Override
		public int compare(Hit_Run r1, Hit_Run r2) {
			String s1 = r1.getReadID();
			String s2 = r2.getReadID();
			return s1.compareTo(s2);
		}

	}

}
