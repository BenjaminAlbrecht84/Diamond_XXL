package util;

public class AlignmentEvaluater {

	public int[] run(String[] alignment, ScoringMatrix matrix) {

		int matches = 0, mismatches = 0, gapOpen = 0, gaps = 0, positives = 0, queryLength = 0;
		String s1 = alignment[0];
		String s2 = alignment[1];
		for (int i = 0; i < s1.length(); i++) {
			char c1 = s1.charAt(i);
			char c2 = s2.charAt(i);
			if (matrix.getScore(c1, c2) > 0)
				positives++;
			if (c1 == c2) {
				matches++;
			} else if (c1 != '-' && c2 != '-') {
				mismatches++;
			} else {
				gaps++;
				if (c1 == '-' && (i == 0 || s1.charAt(i - 1) != '-'))
					gapOpen++;
				else if (c2 == '-' && (i == 0 || s2.charAt(i - 1) != '-'))
					gapOpen++;
			}
			if (c1 != '-')
				queryLength++;
		}

		int[] stats = { alignment[0].length(), matches, positives, mismatches, gapOpen, gaps, queryLength };
		return stats;

	}

}
