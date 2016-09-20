package pipeline.diamond;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import startUp.Main.ParseMode;

public class Diamond_Runner {

	// input options
	private File bin;
	private String optionString;

	// output files
	private File samFile, daaFile;

	// parsed options
	private double lambda, k;
	private int gapOpen, gapExtend;
	private String matrixType;

	public Diamond_Runner(File bin, String optionString) {
		this.bin = bin;
		this.optionString = optionString;
	}

	public File execute(File db, File fna, File out, int cores, ParseMode parseMode) {

		System.out.println("STEP_2>Blasting shredded reads with DIAMOND...");
		long time = System.currentTimeMillis();

		if (blastx(db, fna, out, cores) != 0)
			System.exit(1);

		daaFile = new File(out.getAbsoluteFile() + File.separator + fna.getName().split("\\.")[0] + ".daa");

		long runtime = (System.currentTimeMillis() - time) / 1000;
		System.out.println("OUTPUT>Blasting of the shredded reads with DIAMOND finished. [" + runtime + "s]\n");

		if (parseMode == ParseMode.DAA)
			return daaFile;

		if (view_samFile(daaFile, out) != 0)
			System.exit(1);

		return samFile;

	}

	private int blastx(File db, File fna, File out, int cores) {
		out.mkdir();
		if (db.exists() && fna.exists() && out.exists()) {

			// default setting
			String cmd = bin.getAbsolutePath() + " blastx -d " + db.getAbsolutePath() + " -q " + fna.getAbsolutePath() + " -a "
					+ out.getAbsolutePath() + File.separator + fna.getName().split("\\.")[0] + " " + optionString;

			// sensitive setting
			// String cmd = bin.getAbsolutePath() + " blastx -d " + db.getAbsolutePath() + " -q " + fna.getAbsolutePath() + " -a "
			// + out.getAbsolutePath() + File.separator + fna.getName().split("\\.")[0] + " -p " + cores + " -k 100 --seg no --more-sensitive";

			return executingCommand(cmd);
		} else if (!db.exists())
			System.out.println("Error running Diamond: " + db.getAbsolutePath() + " does not exist!");
		else if (!fna.exists())
			System.out.println("Error running Diamond: " + fna.getAbsolutePath() + " does not exist!");
		else if (out.exists())
			System.out.println("Error running Diamond: " + out.getAbsolutePath() + " does not exist!");

		return 1;
	}

	private int view_samFile(File daaFile, File out) {
		out.mkdir();
		if (daaFile.exists() && out.exists()) {
			samFile = new File(out.getAbsolutePath() + File.separator + daaFile.getName().split("\\.")[0] + ".sam");
			String cmd = bin.getAbsolutePath() + " view -a " + daaFile.getAbsolutePath() + " -o " + out.getAbsolutePath() + File.separator
					+ daaFile.getName().split("\\.")[0] + ".sam -f sam ";
			return executingCommand(cmd);
		} else if (!samFile.exists())
			System.out.println("Error running Diamond: " + samFile.getAbsolutePath() + " does not exist!");
		else if (!out.exists())
			System.out.println("Error running Diamond: " + out.getAbsolutePath() + " does not exist!");

		return 1;
	}

	private int executingCommand(String command) {
		try {

			System.out.println("OUTPUT>Executing " + command);
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec(command);

			// checking error messages
			StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

			// checking error messages
			StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

			errorGobbler.start();
			outputGobbler.start();
			int exitVal = proc.waitFor();

			return exitVal;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return 1;
	}

	class StreamGobbler extends Thread {
		InputStream is;
		String type;

		StreamGobbler(InputStream is, String type) {
			this.is = is;
			this.type = type;
		}

		public void run() {
			try {
				InputStreamReader isr = new InputStreamReader(is);
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("Scoring parameters: ")) {
						String[] columns = line.split(" ");
						matrixType = columns[2].split("=")[1];
						lambda = Double.parseDouble(columns[3].split("=")[1]);
						k = Double.parseDouble(columns[4].split("=")[1]);
						gapOpen = Integer.parseInt(columns[5].split("=")[1].split("/")[0]);
						gapExtend = Integer.parseInt(columns[5].split("=")[1].split("/")[1].split("\\)")[0]);
					}
					System.out.println(type + ">" + line);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
	}

	public double getLambda() {
		return lambda;
	}

	public double getK() {
		return k;
	}

	public String getMatrixType() {
		return matrixType;
	}

	public int getGapOpen() {
		return gapOpen;
	}

	public int getGapExtend() {
		return gapExtend;
	}

	public File getDaaFile() {
		return daaFile;
	}

}
