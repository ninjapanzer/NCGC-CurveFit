/*
<h2 align=center>PUBLIC DOMAIN NOTICE<h2/>

<h3 align=center>NIH Chemical Genomics Center<h3/>
<h3 align=center>National Human Genome Research Institute<h3/>

<p>This software/database is a "United States Government Work" under the 
terms of the United States Copyright Act.  It was written as part of 
the author's official duties as United States Government employee and 
thus cannot be copyrighted.  This software/database is freely 
available to the public for use. The NIH Chemical Genomics Center 
(NCGC) and the U.S. Government have not placed any restriction on its 
use or reproduction. 

<p>Although all reasonable efforts have been taken to ensure the accuracy 
and reliability of the software and data, the NCGC and the U.S. 
Government do not and cannot warrant the performance or results that 
may be obtained by using this software or data. The NCGC and the U.S. 
Government disclaim all warranties, express or implied, including 
warranties of performance, merchantability or fitness for any 
particular purpose. 

<p>Please cite the authors in any work or product based on this material. 
 */

package gov.nih.ncgc.batch;

public class HtsUtil {
	public static boolean isDecimal(String text) throws Exception {
		if (text.length() < 1)
			return false;

		char c = text.charAt(0);
		if (text.length() == 1 && !Character.isDigit(c))
			return false;

		else if (!(Character.isDigit(c) | c == '-' | c == '.'))
			return false;

		for (int i = 1; i < text.length(); i++) {
			c = text.charAt(i);
			if (!(Character.isDigit(c) | c == '.'))
				return false;
		}

		return true;
	}

	public static double[] lineFit(double[] xs, double ys[]) {
		int no = xs.length;

		boolean[] flags = new boolean[no];

		for (int i = 0; i < no; i++)
			flags[i] = true;

		return lineFit(flags, xs, ys);
	}

	public static double[] lineFit(boolean flags[], double[] xs, double ys[]) {
		int no = xs.length;

		double xm = 0;
		double ym = 0;

		for (int i = 0; i < no; i++)
			if (flags[i]) {
				xm += xs[i];
				ym += ys[i];
			}

		xm /= no;
		ym /= no;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;
		for (int i = 0; i < no; i++)
			if (flags[i]) {
				lxx += (xs[i] - xm) * (xs[i] - xm);
				lyy += (ys[i] - ym) * (ys[i] - ym);
				lxy += (xs[i] - xm) * (ys[i] - ym);
			}

		double b = lxy / lxx;
		double a = ym - b * xm;
		double r = lxy / Math.sqrt(lxx * lyy);

		double[] values = { b, a, r };

		return values;
	}

	public static double parseConcentrationMolar(String text)

	{

		text = text.trim().toLowerCase();

		String value = text.substring(0, text.length() - 2);

		double d = Double.parseDouble(value);

		if (text.endsWith("nm"))

			return d * 1e-9;

		else if (text.endsWith("mm"))

			return d * 1e-3;

		else if (text.endsWith("um"))

			return d * 1e-6;

		else if (text.endsWith("m"))

			return Double.parseDouble(text.substring(0, text.length() - 1));
		else
			return Double.parseDouble(value);

	}

}