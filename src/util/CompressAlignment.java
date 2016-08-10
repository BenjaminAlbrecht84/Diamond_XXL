package util;

public class CompressAlignment {

	public String[] run(String[] alignment) {
		
		String cigar = deriveCIGAR(alignment);
		String seq = deriveSEQ(alignment);
		String md = deriveMD(alignment);

		String[] result = { cigar, seq, md };
		return result;

	}

	private String deriveSEQ(String[] ali) {
		String seq = ali[0].replaceAll("-", "");
		return seq;
	}

	private String deriveMD(String[] ali) {
		StringBuffer md = new StringBuffer("");
		char lastType = '-';
		int num = 0;
		for (int i = 0; i < ali[0].length(); i++) {
			char type = getMDType(ali[0].charAt(i), ali[1].charAt(i));
			if (type == 'M')
				num++;
			else if (type != 'I') {
				if (num != 0) {
					md = md.append(num);
					num = 0;
				}
				if (type == 'X') {
					if (lastType == 'D')
						md.append('0');
					md.append(ali[1].charAt(i));
				}
				if (type == 'D') {
					if (lastType != 'D')
						md.append('^');
					md.append(ali[1].charAt(i));
				}
			}
			if (type != 'I')
				lastType = type;
		}
		if (lastType == 'M')
			md = md.append(num);
		return md.toString();
	}

	private char getMDType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		if (c1 != c2)
			return 'X';
		return 'M';
	}

	private String deriveCIGAR(String[] ali) {
		StringBuffer cigar = new StringBuffer("");
		char lastType = getCIGARType(ali[0].charAt(0), ali[1].charAt(0));
		int num = 1;
		for (int i = 1; i < ali[0].length(); i++) {
			char type = getCIGARType(ali[0].charAt(i), ali[1].charAt(i));
			if (type == lastType)
				num++;
			else {
				cigar = cigar.append(num + "" + lastType);
				num = 1;
			}
			lastType = type;
		}
		cigar = cigar.append(num + "" + lastType);
		return cigar.toString();
	}

	private char getCIGARType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		return 'M';
	}

}
