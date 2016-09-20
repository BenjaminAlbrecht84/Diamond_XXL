package startUp;

import io.OptionParser;
import io.daa.DAA_Reader;
import io.daa.DAA_Writer;

import java.io.File;
import java.math.BigInteger;
import java.util.Locale;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import pipeline.diamond.Diamond_Runner;
import pipeline.init.Shotgun_inParallel;
import pipeline.post.HitRun_Rater;
import pipeline.post.HitRun_Writer;
import pipeline.post.HitToSamConverter;
import pipeline.post.Hit_Run;
import pipeline.post.ReadHits;
import pipeline.post.SAM_Parser_inParallel_Binary;
import pipeline.post.mode_one.Alignment_Completer;
import pipeline.post.mode_one.Alignment_Generator_inParallel;
import util.ScoringMatrix;

public class Main {

	public enum ParseMode {
		DAA, SAM
	}

	public static void main(String[] args) {

		ParseMode parseMode = Main.ParseMode.SAM;

		// ensuring English output-style
		Locale.setDefault(new Locale("en", "US"));

		// parsing options from the command line
		OptionParser parser = new OptionParser();
		if (parser.run(args) != 0)
			System.exit(1);

		// running diamond
		int cores = parser.getNumOfThreads();
		File exeFile = parser.getExeFile();
		File dbFile = parser.getDbFile();
		File queryFile = parser.getQueryFile();
		File outputFolder = parser.getOutputFolder();
		String queryName = queryFile.getName().split("\\.")[0];

		// setting filters
		boolean useFilters = true;
		double maxEValue = parser.getMaxSumProbability();
		int minSumScore = parser.getMinSumScore();
		int minCoverage = parser.getMinCoverage();

		// length/step must be a multiple of 3
		int length = adaptSize(parser.getShredLength());
		int step = cmpStep(length, 100 - parser.getShredOverlap());

		// performing sequence shotgun
		Shotgun_inParallel shotgun = new Shotgun_inParallel();
		File shredded_query = shotgun.run(queryFile, outputFolder, length, step, cores);

		// running diamond
		Diamond_Runner dR = new Diamond_Runner(exeFile, parser.getDiamondOptionString());
		File daaFile = null, samFile = null;
		if (parseMode == ParseMode.SAM) {
			samFile = dR.execute(dbFile, shredded_query, outputFolder, cores, parseMode);
			daaFile = dR.getDaaFile();
		} else {
			daaFile = dR.execute(dbFile, shredded_query, outputFolder, cores, parseMode);
			samFile = new File(outputFolder.getAbsolutePath() + File.separator + queryName + "_tmp.sam");
		}

		// parsing db properties
		DAA_Reader daaReader = new DAA_Reader(daaFile);
		BigInteger dbSize = daaReader.getDAAHeader().getDbLetters();

		// parsing sam file
		daaReader = parseMode == ParseMode.DAA ? daaReader : null;
		ScoringMatrix matrix = new ScoringMatrix(dR.getMatrixType(), dR.getGapOpen(), dR.getGapExtend());
		HitRun_Rater scorer = new HitRun_Rater(dR.getLambda(), dR.getK(), dbSize, daaReader, matrix, shotgun.getReadToOrigLength());
		ConcurrentHashMap<String, ReadHits> readMap = null;
		if (parseMode == ParseMode.SAM)
			readMap = new SAM_Parser_inParallel_Binary().parse_Hits(samFile, matrix, scorer, cores);
		else
			readMap = daaReader.parseAllHits(cores);

		// merging/selecting hits
		Vector<Hit_Run> hitRuns = new Alignment_Generator_inParallel(readMap, scorer, samFile, daaReader, matrix, maxEValue, minSumScore, useFilters,
				step, length).run(cores);

		// initializing sam writer
		File out_sam = new File(outputFolder.getAbsolutePath() + File.separator + queryName + ".sam");
		HitToSamConverter samConverter = new HitToSamConverter(samFile, daaReader, step, dbSize, shotgun.getReadToOrigLength(), matrix, out_sam);

		// initializing daa writer
		File out_daa = new File(outputFolder.getAbsolutePath() + File.separator + queryName + ".daa");
		DAA_Writer daaWriter = parseMode == ParseMode.DAA ? new DAA_Writer(out_daa, samFile, daaReader) : null;

		// initializing run writer
		File out_runs = new File(outputFolder.getAbsolutePath() + File.separator + queryName + ".runs");
		HitRun_Writer runWriter = new HitRun_Writer(out_runs, samFile, daaReader, matrix);

		// completing alignments
		new Alignment_Completer().run(hitRuns, queryFile, dbFile, samFile, daaReader, matrix, dR.getLambda(), dR.getK(), cores, scorer, step,
				samConverter, daaWriter, runWriter, maxEValue, minSumScore, minCoverage, useFilters);

		if (parseMode == ParseMode.DAA)
			samFile.delete();

	}

	private static int cmpStep(int length, int p) {
		int step = Math.round(length / p);
		return adaptSize(step);
	}

	private static int adaptSize(int s) {
		while (s % 3 != 0)
			s++;
		return s;
	}

}
