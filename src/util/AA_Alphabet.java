package util;

import java.util.HashMap;

public class AA_Alphabet {

	private final static String aaString = "ARNDCQEGHILKMFPSTWYVBJZX*/\\";

	private static HashMap<Character, Integer> aaToIndex = new HashMap<Character, Integer>();
	static {
		for (int i = 0; i < aaString.length(); i++)
			aaToIndex.put(aaString.charAt(i), i);
	}

	public static String getAaString() {
		return aaString;
	}

	public static char getCharacter(int i) {
		return aaString.charAt(i);
	}

	public static Integer getIndex(char c) {
		if (aaToIndex.containsKey(c))
			return aaToIndex.get(c);
		return null;
	}

}
