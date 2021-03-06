package util;

import java.util.HashMap;
import java.util.Vector;

public class DAACompressAlignment {

	private HashMap<Character, Integer> aaToIndex = new HashMap<Character, Integer>();;

	public Vector<Byte> run(String[] ali) {

		// String aaString = "ARNDCQEGHILKMFPSTWYVBJZX*";
		String aaString = new AA_Alphabet().getAaString();
		for (int i = 0; i < aaString.length(); i++)
			aaToIndex.put(aaString.charAt(i), i);

		Vector<Byte> editOps = new Vector<Byte>();
		char lastType = '-';
		int num = 0;
		for (int i = 0; i < ali[0].length(); i++) {

			char c1 = ali[0].charAt(i);
			char c2 = ali[1].charAt(i);
			char type = getEditType(c1, c2);

			if (type == 'M' || type == 'I') {
				if (type != lastType && num != 0) {
					editOps.addAll(getEditOperation(lastType, num));
					num = 0;
				}
				num++;
			} else {
				if (num != 0 && lastType != '-') {
					editOps.addAll(getEditOperation(lastType, num));
					num = 0;
				}
				editOps.addAll(getEditOperation(type, c2));
			}
			lastType = type;

		}
		if (lastType == 'M' || lastType == 'I')
			editOps.addAll(getEditOperation(lastType, num));

		return editOps;
	}

	private Vector<Byte> getEditOperation(char type, int total) {
			
		Vector<Byte> opVec = new Vector<Byte>();
		while (total > 0) {
			int num = total > 63 ? 63 : total;
			total -= num;
			byte op = 0;
			if (type == 'I')
				op |= 1 << 6;
			op |= num;
			opVec.add(op);
		}
		return opVec;
	}

	private Vector<Byte> getEditOperation(char type, char c) {

		byte op = 0;
		if (type == 'D')
			op |= 2 << 6;
		else
			op |= 3 << 6;

		op |= aaToIndex.get(c);

		Vector<Byte> opVec = new Vector<Byte>();
		opVec.add(op);
		return opVec;
	}

	private char getEditType(char c1, char c2) {
		if (c1 == '-')
			return 'D';
		if (c2 == '-')
			return 'I';
		if (c1 != c2)
			return 'S';
		return 'M';
	}

}
