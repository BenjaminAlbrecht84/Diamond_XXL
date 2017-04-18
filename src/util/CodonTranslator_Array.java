package util;

import java.util.HashMap;

public class CodonTranslator_Array {

	private char[][][] codons;

	public CodonTranslator_Array() {

		codons = new char[26][26][26];
		for (int a = 0; a < 26; a++) {
			for (int b = 0; b < 26; b++) {
				for (int c = 0; c < 26; c++) {
					codons[a][b][c] = 'X';
				}
			}
		}

		addCodon("GGU", 'G');
		addCodon("GGC", 'G');
		addCodon("GGA", 'G');
		addCodon("GGG", 'G');

		addCodon("GCU", 'A');
		addCodon("GCC", 'A');
		addCodon("GCA", 'A');
		addCodon("GCG", 'A');

		addCodon("GUU", 'V');
		addCodon("GUC", 'V');
		addCodon("GUA", 'V');
		addCodon("GUG", 'V');

		addCodon("UUA", 'L');
		addCodon("UUG", 'L');
		addCodon("CUU", 'L');
		addCodon("CUC", 'L');
		addCodon("CUA", 'L');
		addCodon("CUG", 'L');

		addCodon("AUU", 'I');
		addCodon("AUC", 'I');
		addCodon("AUA", 'I');

		addCodon("UCU", 'S');
		addCodon("UCC", 'S');
		addCodon("UCA", 'S');
		addCodon("UCG", 'S');
		addCodon("AGU", 'S');
		addCodon("AGC", 'S');

		addCodon("ACU", 'T');
		addCodon("ACC", 'T');
		addCodon("ACA", 'T');
		addCodon("ACG", 'T');

		addCodon("GAU", 'D');
		addCodon("GAC", 'D');

		addCodon("GAA", 'E');
		addCodon("GAG", 'E');

		addCodon("AAU", 'N');
		addCodon("AAC", 'N');

		addCodon("CAA", 'Q');
		addCodon("CAG", 'Q');

		addCodon("AAG", 'K');
		addCodon("AAA", 'K');

		addCodon("CGU", 'R');
		addCodon("CGC", 'R');
		addCodon("CGA", 'R');
		addCodon("CGG", 'R');
		addCodon("AGA", 'R');
		addCodon("AGG", 'R');

		addCodon("CAU", 'H');
		addCodon("CAC", 'H');

		addCodon("UUU", 'F');
		addCodon("UUC", 'F');

		addCodon("UGU", 'C');
		addCodon("UGC", 'C');

		addCodon("UGG", 'W');

		addCodon("UAU", 'Y');
		addCodon("UAC", 'Y');

		addCodon("AUG", 'M');

		addCodon("CCU", 'P');
		addCodon("CCC", 'P');
		addCodon("CCA", 'P');
		addCodon("CCG", 'P');
	}

	private void addCodon(String codon, char r) {
		int a = getAccess(codon.charAt(0));
		int b = getAccess(codon.charAt(1));
		int c = getAccess(codon.charAt(2));
		codons[a][b][c] = r;
	}

	private int getAccess(char c) {
		return c > 122 ? c - 32 - 65 : c - 65;
	}

	public char translateCodon(String codon) {
		if (codon.length() != 3)
			return 'X';

		int a = getAccess(codon.charAt(0));
		int b = getAccess(codon.charAt(1));
		int c = getAccess(codon.charAt(2));
		return codons[a][b][c];

	}

	public String translate(String dna) {
		String rna = dna.replaceAll("T", "U");
		StringBuffer prot = new StringBuffer("");
		String codon = "";
		for (int i = 0; i < rna.length(); i++) {
			codon = codon.concat("" + rna.charAt(i));
			if (codon.length() == 3) {
				char p = translateCodon(codon);
				prot = prot.append(p);
				codon = "";
			}
		}
		return prot.toString();
	}

	public char translateCodon(char c1, char c2, char c3) {
		// System.out.println(c1 + " " + c2 + " " + c3);
		int a = getAccess(c1);
		int b = getAccess(c2);
		int c = getAccess(c3);
		return codons[a][b][c];
	}

}
