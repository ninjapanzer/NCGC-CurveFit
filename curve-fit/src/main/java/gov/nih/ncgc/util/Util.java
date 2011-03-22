package gov.nih.ncgc.util;

/*

 parse sw.xml and put them into oracle database
 */

import java.text.DecimalFormat;

public class Util {

	private static String[] patches;

	static {
		patches = new String[10];

		for (int i = 1; i < 10; i++) {
			patches[i] = new String();
			for (int j = 0; j < i; j++)
				patches[i] = patches[i] + " ";
		}
	}

	public static String formatInteger(int k, int n) {
		String str = String.valueOf(k);

		return str + patches[n - str.length()];
	}

	public static String formatDouble(double d, int left, int right) {
		String dString = String.valueOf(d);

		int pos = dString.indexOf('.');

		if (pos < 0)
			return dString;

		String leftString = dString.substring(0, pos);

		int delta = dString.length() - pos;

		String rightString = "";

		if (delta > right)
			rightString = dString.substring(pos + 1, pos + 1 + right);
		else
			rightString = dString.substring(pos + 1);

		rightString = "." + rightString;

		delta = left - leftString.length();

		String patch = "";
		for (int i = 0; i < delta; i++)
			patch += " ";

		return patch + leftString + rightString;
	}

	public static String formatDouble(String ds, String pattern) {
		double d = Double.parseDouble(ds);
		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		return output;
	}

	public static String formatDouble(double d, String pattern) {
		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		return output;
	}

	public static String formatDouble(double d) {
		String pattern = "0.#####E0";

		DecimalFormat myFormatter = new DecimalFormat(pattern);

		return myFormatter.format(d);
	}

	public static String formatDouble(double d, int right) {

		String pattern = "#0.";
		for (int i = 0; i < right; i++)
			pattern = pattern + "0";

		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		if (true)
			return output;

		// DecimalFormat format = new DecimalFormat("#,##0.00000000");
		// FieldPosition f = new FieldPosition(0);
		// StringBuffer s = new StringBuffer();
		// String dString = format.format(Double.parseDouble(d), s,
		// f).toString();

		String dString = String.valueOf(d);

		int pos = dString.indexOf('.');

		if (pos < 0)
			return dString;

		String leftString = dString.substring(0, pos);

		int delta = dString.length() - pos;

		String rightString = "";

		if (delta > right)
			rightString = dString.substring(pos + 1, pos + 1 + right);
		else
			rightString = dString.substring(pos + 1);

		rightString = "." + rightString;

		return leftString + rightString;
	}

	public static String formatSeq(String q0, String q1, String s0, String s1, String qSeq, String sSeq, String mSeq) {
		int qFrom = Integer.parseInt(q0);
		int qTo = Integer.parseInt(q1);

		int sFrom = Integer.parseInt(s0);
		int sTo = Integer.parseInt(s1);

		StringBuffer buffer = new StringBuffer();

		int start = 0;
		int end = 0;

		int len = qSeq.length() / 60 + 1;
		for (int i = 0; i < len; i++) {
			int k = i * 60;
			int l = (i + 1) * 60;
			if (l > qSeq.length())
				l = qSeq.length();

			// quer
			if (qFrom < qTo) {
				start = qFrom + k;
				end = qFrom + l;
			} else {
				start = qFrom - k;
				end = qFrom - l;
			}

			buffer.append("Query:   ");
			buffer.append(formatInteger(start, 10));
			buffer.append(qSeq.substring(k, l));
			buffer.append(" ");
			buffer.append(formatInteger(end, 10));
			buffer.append("\n");

			// middle line
			buffer.append("         ");
			buffer.append("          ");
			buffer.append(mSeq.substring(k, l));
			buffer.append("\n");

			// subject
			if (sFrom < sTo) {
				start = sFrom + k;
				end = sFrom + l;
			} else {
				start = sFrom - k;
				end = sFrom - l;
			}

			buffer.append("Subject: ");
			buffer.append(formatInteger(start, 10));
			buffer.append(sSeq.substring(k, l));
			buffer.append(" ");
			buffer.append(formatInteger(end, 10));
			buffer.append("\n");
			buffer.append("\n");
		}

		return buffer.toString();
	}

}
