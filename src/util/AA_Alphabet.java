package util;

public class AA_Alphabet {

	private final static String aaString = "ARNDCQEGHILKMFPSTWYVBJZX*";

	public String getAaString() {
		return aaString;
	}

	public char getCharacter(int i) {
		return aaString.charAt(i);
	}

}
