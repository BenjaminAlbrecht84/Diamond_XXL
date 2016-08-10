package io.daa;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class DAA_Header {

	private File daaFile;

	// header one
	protected long magicNumber;
	protected long version;

	// header two
	protected BigInteger dbLetters;
	protected long diamondBuild, dbSeqs, dbSeqsUsed, flags, queryRecords;
	protected int modeRank, gapOpen, gapExtend, reward, penalty, reserved1, reserved2, reserved3;

	protected double k, lambda, reserved4, reserved5;
	protected final byte[] scoreMatrix = new byte[16];
	protected final long[] blockSize = new long[256];
	protected final byte[] blockTypeRank = new byte[256];

	// reference information
	private byte[][] references;
	protected int[] refLengths;

	private final int referenceLocationChunkBits = 6; // 6 bits = 64 chunk size
	private final int referenceLocationChunkSize = 1 << referenceLocationChunkBits;
	private long[] referenceLocations; // location of every
										// 2^referenceLocationChunkBits
										// reference

	// reference annotations:
	protected int numberOfRefAnnotations;
	protected int[][] refAnnotations = new int[256][];
	protected String[] refAnnotationNames = new String[256];
	protected int refAnnotationIndexForTaxonomy = -1;

	// helper variables:
	protected String scoreMatrixName;

	protected long headerSize;

	protected int refNamesBlockIndex = -1;
	protected int refLengthsBlockIndex = -1;
	protected int alignmentsBlockIndex = -1;

	public DAA_Header(File daaFile) {
		this.daaFile = daaFile;
		if (magicNumber == 0)
			load();
	}

	private void load() {

		try {

			InputStream is = new BufferedInputStream(new FileInputStream(daaFile));
			ByteBuffer buffer = ByteBuffer.allocate(2448);
			buffer.order(ByteOrder.LITTLE_ENDIAN);
			is.read(buffer.array());

			try {

				// magicNumber = buffer.getLong(0);
				// version = buffer.getLong(8);

				// diamondBuild = buffer.getLong(16);
				// dbSeqs = buffer.getLong(24);
				dbSeqsUsed = buffer.getLong(32);
				dbLetters = BigInteger.valueOf(buffer.getLong(40));
				// flags = buffer.getLong(48);
				queryRecords = buffer.getLong(56);

				modeRank = buffer.getInt(64);
				gapOpen = buffer.getInt(68);
				gapExtend = buffer.getInt(72);
				reward = buffer.getInt(76);
				penalty = buffer.getInt(80);
				reserved1 = buffer.getInt(84);
				reserved2 = buffer.getInt(88);
				reserved3 = buffer.getInt(92);

				k = buffer.getDouble(96);
				lambda = buffer.getDouble(104);
				reserved4 = buffer.getDouble(112);
				reserved5 = buffer.getDouble(120);

				// buffer.position(128);
				// for (int i = 0; i < scoreMatrix.length; i++)
				// scoreMatrix[i] = (byte) buffer.get();

				buffer.position(144);
				for (int i = 0; i < blockSize.length; i++)
					blockSize[i] = buffer.getLong();

				buffer.position(2192);
				for (int i = 0; i < blockTypeRank.length; i++) {
					blockTypeRank[i] = (byte) buffer.get();
					switch (blockTypeRank[i]) {
					case 1:
						alignmentsBlockIndex = i;
						break;
					case 2:
						refNamesBlockIndex = i;
						break;
					case 3:
						refLengthsBlockIndex = i;
						break;
					}
				}

				buffer.position(2448);
				headerSize = buffer.position();

			} finally {
				is.close();
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void loadAllReferences() {
		try {
			RandomAccessFile raf = new RandomAccessFile(daaFile, "r");
			try {

				raf.seek(getLocationOfBlockInFile(refNamesBlockIndex));

				references = new byte[(int) dbSeqsUsed][];
				referenceLocations = new long[1 + ((int) dbSeqsUsed >>> referenceLocationChunkBits)];
				for (int i = 0; i < (int) dbSeqsUsed; i++) {
					if ((i & (referenceLocationChunkSize - 1)) == 0)
						referenceLocations[i >>> referenceLocationChunkBits] = raf.getFilePointer();
					int c = raf.read();
					while (c != 0)
						c = raf.read();
				}

				refLengths = new int[(int) dbSeqsUsed];
				for (int i = 0; i < dbSeqsUsed; i++) {
					ByteBuffer buffer = ByteBuffer.allocate(4);
					buffer.order(ByteOrder.LITTLE_ENDIAN);
					raf.read(buffer.array());
					refLengths[i] = buffer.getInt();
				}

			} finally {
				raf.close();
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public byte[] getReference(int index) throws IOException {
		if (references[index] != null)
			return references[index];

		RandomAccessFile raf = new RandomAccessFile(daaFile, "r");

		try {

			int iChunk = (index >>> referenceLocationChunkBits);
			raf.seek(referenceLocations[iChunk]);

			int start = iChunk * referenceLocationChunkSize;
			int stop = Math.min((int) dbSeqsUsed, start + referenceLocationChunkSize);

			for (int i = start; i < stop; i++) {
				StringBuffer buf = new StringBuffer();
				byte b = (byte) raf.read();
				while (b != 0) {
					buf.append((char) b);
					b = (byte) raf.read();
				}
				references[i] = buf.toString().getBytes();
			}

		} finally {
			raf.close();
		}

		return references[index];
	}

	public int getAlignmentsBlockIndex() {
		return alignmentsBlockIndex;
	}

	public long getLocationOfBlockInFile(int blockIndex) {
		long location = headerSize;
		for (int i = 0; i < blockIndex; i++)
			location += blockSize[i];
		return location;
	}

	public long getNumberOfQueryRecords() {
		return queryRecords;
	}

	public int getRefLength(int i) {
		return refLengths[i];
	}

	public double getK() {
		return k;
	}

	public double getLambda() {
		return lambda;
	}

	public BigInteger getDbLetters() {
		return dbLetters;
	}

}
