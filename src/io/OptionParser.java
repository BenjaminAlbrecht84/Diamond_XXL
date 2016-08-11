package io;

import java.io.File;

public class OptionParser {

	private int numOfThreads = Runtime.getRuntime().availableProcessors(), length = 1000;
	private File dbFile = null, queryFile = null, alignmentFile = null, exeFile = null;
	private double maxEValue = 0.001;
	private int minSumScore = 30;

	public int run(String[] args) {

		if (args.length == 1 && (args[0] == "help" || args[0] == "-h")) {
			printOptions();
			System.exit(0);
		}

		for (int i = 0; i < args.length - 1; i++) {
			String option = args[i];
			switch (option) {
			case "-s":
				try {
					minSumScore = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "-e":
				try {
					maxEValue = Double.parseDouble(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a double " + (args[i + 1]));
				}
				i++;
				break;
			case "-p":
				try {
					int p = Integer.parseInt(args[i + 1]);
					numOfThreads = p < numOfThreads ? p : numOfThreads;
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "-exe":
				try {
					exeFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "-d":
				try {
					dbFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "-q":
				try {
					queryFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "-a":
				try {
					alignmentFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "-l":
				try {
					length = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			default:
				System.out.println("ERROR: Uknown option " + option);
			}
		}

		if (exeFile == null)
			System.err.println("ERROR: missing path to the DIAMOND executable (Option -exe)");
		if (dbFile == null)
			System.err.println("ERROR: missing path DIAMOND database file (Option -d)");
		if (queryFile == null)
			System.err.println("ERROR: missing query input file in FASTA or FASTQ format (Option -q)");
		if (alignmentFile == null)
			System.err.println("ERROR: missing path to output file in SAM format (Option -a)");

		if (exeFile == null || dbFile == null || queryFile == null || alignmentFile == null || length == -1) {
			printOptions();
			return 1;
		}

		if (!exeFile.exists())
			System.err.println("ERROR: the DIAMOND executable does not exist " + exeFile.getAbsolutePath());
		if (!dbFile.exists())
			System.err.println("ERROR: DIAMOND database file does not exist " + dbFile.getAbsolutePath());
		if (!queryFile.exists())
			System.err.println("ERROR: query input file does not exist " + queryFile.getAbsolutePath());
		if (!alignmentFile.exists() && !alignmentFile.mkdir())
			System.err.println("ERROR: not a valid path for the output file " + alignmentFile.getAbsolutePath());

		if (!exeFile.exists() || !dbFile.exists() || !queryFile.exists() || (!alignmentFile.exists() && !alignmentFile.mkdir())) {
			printOptions();
			return 1;
		}

		return 0;

	}

	private void printOptions() {
		System.out.println("Mandatory Options:");
		System.out.println("-exe \t Path to the DIAMOND executable.");
		System.out.println("-d \t Path to DIAMOND database file (not including the file extension .dmnd).");
		System.out.println("-q \t Path to query input file in FASTA or FASTQ format (may be gzip compressed).");
		System.out.println("-a \t Path to output file in SAM format (extension .daa will be appended).");
	}

	public File getExeFile() {
		return exeFile;
	}

	public int getNumOfThreads() {
		return numOfThreads;
	}

	public File getDbFile() {
		return dbFile;
	}

	public File getQueryFile() {
		return queryFile;
	}

	public File getAlignmentFile() {
		return alignmentFile;
	}

	public int getLength() {
		return length;
	}

	public double getMaxEValue() {
		return maxEValue;
	}

	public int getMinSumScore() {
		return minSumScore;
	}

}
