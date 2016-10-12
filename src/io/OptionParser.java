package io;

import java.io.File;

public class OptionParser {

	// mandatory options
	private File dbFile = null, queryFile = null, outputFolder = null, exeFile = null;

	// performance options
	private int numOfThreads = Runtime.getRuntime().availableProcessors();
	private int shredLength = 1000;
	private int shredOverlap = 90;
	private boolean realign = false;

	// scoring & reporting options
	private Integer delta = null;
	private double maxSumProbability = 0.001;
	private int minSumScore = 30;
	private int minCoverage = 90;

	// diamond options
	private DiamondOptionsContainer diamondOpts = new DiamondOptionsContainer();

	public int run(String[] args) {

		if (args.length == 1 && (args[0] == "help" || args[0] == "-h")) {
			printOptions();
			System.exit(0);
		}

		for (int i = 0; i < args.length; i++) {
			String option = args[i];
			switch (option) {
			case "--minSumScore":
				try {
					minSumScore = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "-maxSumProbability":
				try {
					maxSumProbability = Double.parseDouble(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a double " + (args[i + 1]));
				}
				i++;
				break;
			case "--threads":
			case "-p":
				try {
					int p = Integer.parseInt(args[i + 1]);
					diamondOpts.setThreads(p);
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
			case "--db":
			case "-d":
				try {
					dbFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "--query":
			case "-q":
				try {
					queryFile = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "--out":
			case "-o":
				try {
					outputFolder = new File(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not a file " + (args[i + 1]));
				}
				i++;
				break;
			case "--shred-length":
			case "-l":
				try {
					shredLength = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--shred-overlap":
				try {
					shredOverlap = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--delta":
				try {
					delta = Integer.parseInt(args[i + 1]);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--realign":
				realign = true;
				break;
			case "--salltitles":
				diamondOpts.setSalltitles(true);
				break;
			case "--sensitive":
				diamondOpts.setSensitivityDegree(1);
				break;
			case "--more-sensitive":
				diamondOpts.setSensitivityDegree(2);
				break;
			case "--gapopen":
				try {
					int gapOpen = Integer.parseInt(args[i + 1]);
					diamondOpts.setGapOpen(gapOpen);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--gapextend":
				try {
					int gapextend = Integer.parseInt(args[i + 1]);
					diamondOpts.setGapExtend(gapextend);
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--matrix":
				String matrix = args[i + 1];
				boolean isValid = diamondOpts.setMatrix(matrix);
				if (!isValid)
					System.err.print("ERROR: not a valid scoring matrix " + (args[i + 1]));
				i++;
				break;
			case "--seg":
				String s = args[i + 1];
				isValid = diamondOpts.setSeg(s);
				if (!isValid)
					System.err.print("ERROR: not a valid seq option (yes/no) " + (args[i + 1]));
				i++;
				break;
			case "--maxtarget":
				try {
					diamondOpts.setMaxTargetSeqs(Integer.parseInt(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--top":
				diamondOpts.setTop(true);
				break;
			case "--evalue":
			case "-e":
				try {
					diamondOpts.seteValue(Double.parseDouble(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not a double " + (args[i + 1]));
				}
				i++;
				break;
			case "--minscore":
				try {
					diamondOpts.setMinScore(Integer.parseInt(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--queryCover":
				try {
					diamondOpts.setQueryCover(Integer.parseInt(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--block-size":
				try {
					diamondOpts.setBlockSize(Double.parseDouble(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not a double " + (args[i + 1]));
				}
				i++;
				break;
			case "--index-chunks":
				try {
					diamondOpts.setIndexChunks(Integer.parseInt(args[i + 1]));
				} catch (Exception e) {
					System.err.print("ERROR: not an integer " + (args[i + 1]));
				}
				i++;
				break;
			case "--tmpdir":
				isValid = diamondOpts.setTmpDir(args[i + 1]);
				if (!isValid)
					System.err.print("ERROR: not a valid directory " + (args[i + 1]));
				i++;
				break;
			default:
				System.out.println("ERROR: Uknown option " + option);
				System.exit(1);
			}
		}

		if (exeFile == null)
			System.err.println("ERROR: missing path to the DIAMOND executable (Option -exe)");
		if (dbFile == null)
			System.err.println("ERROR: missing path DIAMOND database file (Option --db/-d)");
		if (queryFile == null)
			System.err.println("ERROR: missing query input file in FASTA or FASTQ format (Option --query/-q)");
		if (outputFolder == null)
			System.err.println("ERROR: missing path to output folder (Option --out/-o)");

		if (exeFile == null || dbFile == null || queryFile == null || outputFolder == null || shredLength == -1) {
			printOptions();
			return 1;
		}

		if (!exeFile.exists())
			System.err.println("ERROR: the DIAMOND executable does not exist " + exeFile.getAbsolutePath());
		if (!dbFile.exists())
			System.err.println("ERROR: DIAMOND database file does not exist " + dbFile.getAbsolutePath());
		if (!queryFile.exists())
			System.err.println("ERROR: query input file does not exist " + queryFile.getAbsolutePath());
		if (!outputFolder.exists() && !outputFolder.mkdir())
			System.err.println("ERROR: not a valid path for the output folder " + outputFolder.getAbsolutePath());

		if (!exeFile.exists() || !dbFile.exists() || !queryFile.exists() || (!outputFolder.exists() && !outputFolder.mkdir())) {
			printOptions();
			return 1;
		}

		return 0;

	}

	private void printOptions() {
		System.out.println("Mandatory Options:");
		System.out.println("-exe \t Path to the DIAMOND binary file.");
		System.out.println("--db/-d \t Path to DIAMOND database file.");
		System.out.println("--query/-q \t Path to query input file in FASTA or FASTQ format.");
		System.out.println("--out/-o \t Path to output folder.");
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

	public File getOutputFolder() {
		return outputFolder;
	}

	public int getShredLength() {
		return shredLength;
	}

	public double getMaxSumProbability() {
		return maxSumProbability;
	}

	public int getMinSumScore() {
		return minSumScore;
	}

	public int getMinCoverage() {
		return minCoverage;
	}

	public Integer getDelta() {
		return delta;
	}

	public int getShredOverlap() {
		return shredOverlap;
	}

	public String getDiamondOptionString() {
		return diamondOpts.getDiamondOptionsString();
	}

	public boolean doRealign() {
		return realign;
	}

}
