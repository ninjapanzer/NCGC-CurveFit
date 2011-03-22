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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math.stat.regression.SimpleRegression;

public class HillFit extends HillConstants {

	public static double LN10 = Math.log(10.0);

	public static boolean debug = false;

	public static boolean fastFlag = false;

	public static boolean info = false;

	private static boolean[] maskings_ = null;

	public static void main(String[] args) {
		try {
			double max = 2510548;

			List list = FileUtils.readLines(new File("c:\\temp\\biphase.txt"));
//			ArrayList list = FileUtil.fileToArrayList("c:\\temp\\biphase.txt");

			double[] xs = new double[list.size()];
			double[] ys = new double[list.size()];

			for (int i = 0; i < list.size(); i++) {
				String line = (String) list.get(i);
				int k = line.indexOf("\t");
				if (k < 0)
					k = line.indexOf(" ");

				xs[i] = Math.log10(Double.parseDouble(line.substring(0, k).trim()));
				ys[i] = 100.0 * Double.parseDouble(line.substring(k).trim()) / max;
			}

			double[] values = HillFit.hillFitFast(HillFit.PI_MAX, HillFit.PS_MIN, null, null, xs, ys, "biphase");

			for (int i = 0; i < values.length; i++)
				System.out.println(values[i]);
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}

	public HillFit() {
		super();
	}

	public static void doHill(String[] args) throws Exception {

		List rowList = FileUtils.readLines(new File(args[0]));
//		ArrayList rowList = FileUtil.fileToArrayList(args[0]);
		int start = Integer.parseInt(args[1]);
		int end = Integer.parseInt(args[2]);
		java.io.PrintStream printStream = System.out;
		if (args.length > 3)
			printStream = new java.io.PrintStream(new java.io.File(args[3]));

		String line = (String) rowList.get(0);

		String[] labels = line.split("\t");
		printStream.println("curve class\t" + line);

		CurveClass cc = new CurveClass();
		cc.setAllowBell(true);

		for (int i = 1; i < rowList.size(); i++) {
			line = (String) rowList.get(i);

			String[] values = line.split("\t");

			HillCurve hc = doHillCurve(start, end, labels, values);

			if (!hc.hasCurve()) {
				printStream.print(cc.enzCurveClass(hc) + "\t");
				printStream.println(line);
				continue;
			}

			double[] fitValues = hc.getObjs();
			double y0 = fitValues[6];
			double x05 = fitValues[0];
			double yinf = fitValues[1];
			double slope = fitValues[2];
			double r2 = fitValues[3];

			printStream.print(cc.enzCurveClass(hc) + "\t");
			printStream.print(line);
			printStream.print("\t" + x05 + "\t" + slope + "\t" + y0 + "\t" + yinf + "\t" + r2);
			if (maskings_ != null)
				for (int j = 0; j < maskings_.length; j++)
					printStream.print("\t" + (maskings_[j] == true ? 1 : 0));
			printStream.println();

		}
	}

	public static double[] doHill(int start, int end, String[] labels, String[] values) throws Exception {
		return doHillCurve(start, end, labels, values).getObjs();
	}

	public static HillCurve doHillCurve(int start, int end, String[] labels, String[] values) throws Exception {
		ArrayList listx = new ArrayList();
		ArrayList listy = new ArrayList();

		for (int i = start; i < labels.length && i < end; i++) {
			if (labels[i] == null || labels[i].equals("null"))
				throw new Exception("Missing labels");

			if (values[i] != null && HtsUtil.isDecimal(values[i])) {
				double conc = HtsUtil.parseConcentrationMolar(labels[i]);

				listx.add(String.valueOf(Math.log10(conc)));
				;

				listy.add(values[i]);
			}
		}

		double[] x = new double[listx.size()];
		double[] y = new double[listy.size()];

		for (int i = 0; i < listx.size(); i++) {
			x[i] = Double.parseDouble((String) listx.get(i));
			y[i] = Double.parseDouble((String) listy.get(i));

			// System.out.println(x[i]+":"+y[i]);
		}

		if (true) {
			Object[] objs = doHill(x, y);

			maskings_ = null;
			if (objs == null)
				return new HillCurve(y, x);
			// return null;

			maskings_ = (boolean[]) objs[0];
			HillCurve hc = new HillCurve(y, x, maskings_, values[0]);
			double[] fitValues = (double[]) objs[1];
			int iterNo = Integer.valueOf((String) objs[3]);
			hc.setCurve(fitValues[0], fitValues[2], fitValues[6], fitValues[1], fitValues[3], iterNo);
			hc.setObjs(fitValues);
			return hc;
		}

		boolean[] flags = maskDif(x, y, null, null, null);

		double[] fitValues = hillFit(PI_MAX, PS_MIN, flags, null, x, y);

		if (fitValues == null) {
			return null;
		}

		double x05 = fitValues[0];
		double yinf = fitValues[1];
		double slope = fitValues[2];
		double r2 = fitValues[3];

		int mask = 0;
		for (int i = 0; i < x.length; i++) {
			if (flags[i] == false)
				mask++;
		}

		// System.out.println(values[0]+"\t"+x05+"\t"+slope+"\t"+yinf+"\t"+r2+"\t"+mask);

		// return fitValues;

		HillCurve hc = new HillCurve(y, x);
		hc.setCurve(fitValues[0], fitValues[2], fitValues[6], fitValues[1], fitValues[3], 0);
		hc.setObjs(fitValues);

		return hc;
	}

	public static boolean[] parseMaskFlag(String maskFlag) {
		if (maskFlag == null)
			return null;

		char[] cols = maskFlag.toCharArray();

		boolean[] flags = new boolean[cols.length];

		for (int i = 0; i < cols.length; i++) {
			if (cols[i] == '1')
				flags[i] = false;
			else
				flags[i] = true; // true means not masked
		}

		return flags;
	}

	public static Object[] doHill(double[] x, double[] y, String maskFlag, int iterNo, String option) {
		boolean[] maskFlags = parseMaskFlag(maskFlag);

		return doHill(x, y, maskFlags, iterNo, option);
	}

	public static Object[] doHill(double[] x, double[] y) {
		return doHill(x, y, HillConstants.PARTIAL_FIT_MASK_NO, true);
	}

	public static Object[] doHill(double[] x, double[] y, int iterNo, boolean p4Fit) {
		return doHill(x, y, null, iterNo, p4Fit);
	}

	public static Object[] doHill(double[] x, double[] y, String maskFlag, int iterNo, boolean p4Fit) {
		return doHill(x, y, maskFlag, iterNo, String.valueOf(p4Fit));
	}

	public static Object[] doHill(double[] x1, double[] y1, boolean[] maskFlags, int iterNo, String option) {
		return doHill(x1, y1, maskFlags, iterNo, option, true);
	}

	public static Object[] doHill(double[] x1, double[] y1, String maskFlag, int iterNo, String option, boolean pctFlag) {
		boolean[] maskFlags = parseMaskFlag(maskFlag);

		return doHill(x1, y1, maskFlags, iterNo, option, pctFlag);
	}

	public static Object[] doHill(double[] x1, double[] y1, boolean[] maskFlags, int iterNo, String option,
			boolean pctFlag) {
		if (option == null)
			option = "";

		if (option.equals("constant_fit"))
			return null;

		// if (HillConstants.bellMask == false && isOnePoint(y1))
		// return null;

		if (x1.length < 3)
			return null;

		int no = x1.length;

		double[] x = new double[no];
		double[] y = new double[no];

		if (x1[0] > x1[no - 1]) {
			for (int i = 0; i < no; i++) {
				x[i] = x1[no - i - 1];
				y[i] = y1[no - i - 1];
			}
		} else {
			for (int i = 0; i < no; i++) {
				x[i] = x1[i];
				y[i] = y1[i];
			}
		}

		// int iterNo = x.length-HillConstants.PARTIAL_FIT_NO+1;

		double yCoef = 1.0;
		// if not pct data, convert into pct
		if (pctFlag == false) {
			double max = -Double.MAX_VALUE;
			for (int i = 0; i < no; i++)
				max = Math.max(max, y[i]);

			yCoef = 0.01 * max;

			for (int i = 0; i < no; i++)
				y[i] /= yCoef;
		}

		double bestR2 = 0.0;
		HashMap[] bestMaps = null;
		double[] bestFitValues = null;
		boolean[] bestFlags = null;
		double bestDelta = 0.0;
		int bestIterNo = 0;

		if (maskFlags != null) {
			iterNo = 1;
		}

		for (int i = 0; i < iterNo; i++) {
			double[] xs = new double[x.length - i];
			double[] ys = new double[xs.length];

			for (int j = 0; j < xs.length; j++) {
				xs[j] = x[j];
				ys[j] = y[j];
			}

			HashMap[] maps = new HashMap[xs.length];

			boolean noChangeFlag = false;
			if (option.toLowerCase().startsWith("bi"))
				noChangeFlag = true;

			boolean[] flags = null;

			if (maskFlags != null) {
				flags = maskFlags;
				noChangeFlag = true;
			} else
				flags = HillFit.maskDif(xs, ys, null, null, maps);

			if (xs.length < 3)
				break;

			double[] fitValues = null;

			if (option.equals("true") || option.equals("false")) {
				boolean p4Fit = Boolean.valueOf(option);

				fitValues = HillFit.hillFit(PI_MAX, PS_MIN, flags, null, xs, ys, noChangeFlag, p4Fit);
			} else {
				fitValues = HillFit.hillFit(PI_MAX, PS_MIN, flags, null, xs, ys, noChangeFlag, option);
			}

			// in case no fit

			if (fitValues == null && bestFlags == null) {
				bestFlags = flags;
			}

			if (fitValues == null)
				continue;

			double r2 = fitValues[3];
			double delta = Math.abs(fitValues[1] - fitValues[6]);

			if (option.toLowerCase().startsWith("bi")) {
				// {y0_min,yi_min,f_min,e1_min,s1_min,e2_min,s2_min,r2,dev_const,dev_min};
				r2 = fitValues[7];
				delta = Math.abs(fitValues[0] - fitValues[1]);
			}

			if (r2 * xs.length > bestR2 && delta > bestDelta) {
				bestR2 = r2 * xs.length;
				bestFitValues = fitValues;
				bestFlags = flags;
				bestMaps = maps;
				bestDelta = delta;
				bestIterNo = i;

				if (info)
					System.out.println(i + "\t" + xs.length + "\t" + r2 + "\t" + bestR2 + "\t" + delta + "\t"
							+ bestDelta);
			} else {
				if (info)
					System.out.println(i + "\t" + xs.length + "\t" + r2 + "\t" + bestR2 + "\t" + delta + "\t"
							+ bestDelta);
				break;
			}
		}

		if (bestFitValues == null) {
			Object[] values = { bestFlags, bestFitValues, bestMaps, String.valueOf(bestIterNo), String.valueOf(yCoef) };

			return values;
		}

		boolean[] flags = new boolean[x.length];
		for (int i = 0; i < bestFlags.length; i++)
			flags[i] = bestFlags[i];

		for (int i = bestFlags.length; i < x.length; i++)
			flags[i] = false;

		HashMap[] maps = new HashMap[x.length];
		for (int i = 0; i < bestFlags.length; i++)
			maps[i] = bestMaps[i];

		for (int i = bestFlags.length; i < x.length; i++)
			maps[i] = null;

		// {y0_min,yi_min,f_min,e1_min,s1_min,e2_min,s2_min,r2,dev_const,dev_min};
		// double[] values = {x05_min, yinf_min, slope_min, r2, dev_const,
		// dev_min, 0.0};
		if (pctFlag == false) {
			bestFitValues[1] *= yCoef;
			bestFitValues[2] *= yCoef;
			bestFitValues[6] *= yCoef;
		}

		Object[] values = { flags, bestFitValues, maps, String.valueOf(bestIterNo), String.valueOf(yCoef) };

		return values;

	}

	public static boolean[] mask(double[] xs, double[] ys) {
		return maskDif(xs, ys, null, null, null);
	}

	public static boolean[] maskDif2(double[] xs, double[] ys, int[][] ds, double[] ms) {
		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		// calc mean and std
		int no = 0;

		double mean = 0.0;
		for (int i = 0; i < xs.length; i++) {
			mean += ys[i];
			no++;
		}
		mean /= no;

		double dev = 0.0;
		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - mean;
			dev += delta * delta;
		}
		dev = Math.sqrt(dev / no);

		// check first point
		// if (false && P0 >0.0 && Math.abs(ys[0]) > P0)

		if (PMAX > 0.0 && Math.abs(ys[0]) > PMAX) {
			if (info)
				System.out.println("Abs(y0)>Pref.PMAX");
			// System.out.println("Abs(y0)>Pref.P0");

			flags[0] = false;
		}

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++) {
			// check against absolute value
			if (PMAX > 0.0 && Math.abs(ys[i]) > PMAX) {
				if (info)
					System.out.println("Abs(y)>Pref.PMAX");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}

			// check against deviation from mean
			double ratio = Math.abs(ys[i] - mean) / dev;

			if (TR > 0.0 && ratio > TR) {
				if (info)
					System.out.println("Deviation from mean > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}

			// shape like v
			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
				vFlag = true;

			// check against angle and delta-the difference from the
			// extrapolated line
			double t1 = 0.0;
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;
			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;
			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta1 = 10000000.0;
			double delta2 = 10000000.0;

			if (i == no) {
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2])
						/ (xs[no - 1] - xs[no - 2]) - ys[no - 2]);
			} else {
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
						- ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;
			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			if (debug)
				System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no) {
				if (delta > TPN && TPN > 0.0) {
					if (info)
						System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				}

				continue;
			}

			if (i != no && delta > TP && TP > 0.0) {
				if (info)
					System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			} else if (vFlag && t < THETA && delta > TPV && THETA > 0.0 && TPV > 0.0) {
				if (info)
					System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"
							+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
			}

			if (ds != null) {
				if (delta > 499)
					delta = 499;
				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++) {
			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;

			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++) {
				String line = (String) list.get(j);
				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));

				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd) {
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}

		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);

		return flags;
	}

	public static boolean[] maskDif(double[] xs, double[] ys, int[][] ds, double[] ms, HashMap[] maps) {
		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		// one point masking
		double sd = HillConstants.CLASSIFICATION_SD;
		double sdf = HillConstants.CLASSIFICATION_SD_FACTOR;

		int pos = -1;
		int no = 0;

		for (int i = 0; i < ys.length; i++) {
			if (Math.abs(ys[i]) >= sd * sdf) {
				no++;
				pos = i;
			}
		}

		if (no == 1 && pos < ys.length / 2) {
			flags[pos] = false;

			return flags;
		}

		// calc mean and std
		no = 0;

		double mean = 0.0;
		double minY = 1000000.0;
		double maxY = -1000000.0;

		for (int i = 0; i < xs.length; i++) {
			maxY = Math.max(maxY, ys[i]);
			minY = Math.min(minY, ys[i]);

			mean += ys[i];

			no++;
		}
		mean /= no;

		double range = maxY - minY;

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - mean;
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / no);

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		// check first point
		// System.out.println(ys[0]+":"+ys[1]);

		if (Math.abs(ys[0]) > HillConstants.P0 && Math.abs(ys[0]) - Math.abs(ys[1]) > HillConstants.MIN_Y_RANGE) {
			if (info)
				System.out.println("abs(y0) > abs(y1) " + ys[0] + ":" + ys[1]);

			flags[0] = false;

			if (maps != null) {
				HashMap map = new HashMap();

				map.put("sample_conc", String.valueOf(xs[0]));
				map.put("sample_resp", String.valueOf(ys[0]));
				map.put("mask_value", String.valueOf(ys[0]));
				map.put("mask_reason", "abs(y0) > abs(y1):" + String.valueOf(ys[0]) + " " + String.valueOf(ys[1]));
				maps[0] = map;
			}
		}

		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++) {
			// check against absolute value; disabled
			if (false && PMAX > 0.0 && Math.abs(ys[i]) > PMAX) {
				if (info)
					System.out.println("Abs(y)>Pref.PMAX");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(ys[i]));
					map.put("mask_reason", "abs(y)>Pref.PMAX:" + String.valueOf(ys[i]) + " " + String.valueOf(PMAX));

					maps[i] = map;
				}

				continue;
			}

			// check against deviation from mean; make it stringent
			double ratio = Math.abs(ys[i] - mean) / dev;

			if (TR > 0.0 && ratio > TR) {
				if (info)
					System.out.println("(resp-mean)/std > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(ratio));
					map.put("mask_reason", "(resp-mean)/std > Pref.TR:" + String.valueOf(ratio) + " "
							+ String.valueOf(TR));

					maps[i] = map;
				}

				continue;
			}

			// shape like v
			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0
					&& Math.abs(ys[i] - ys[i - 1]) > HillConstants.V_DEPTH
					&& Math.abs(ys[i] - ys[i + 1]) > HillConstants.V_DEPTH)
				vFlag = true;

			// check against angle and delta-the difference from the
			// extrapolated line
			double t1 = 0.0;
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;
			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;
			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta1 = 10000000.0;
			double delta2 = 10000000.0;

			if (i == no) {
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2])
						/ (xs[no - 1] - xs[no - 2]) - ys[no - 2]);
			} else {
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
						- ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;
			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			if (debug)
				System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no) {
				if (delta * 100.0 / range > TPN && TPN > 0.0) {
					if (info)
						System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

					if (maps != null) {
						HashMap map = new HashMap();

						map.put("sample_conc", String.valueOf(xs[i]));
						map.put("sample_resp", String.valueOf(ys[i]));
						map.put("mask_value", String.valueOf(delta));
						map.put("mask_reason", "Last point's delta distance > Pref.TPN:" + String.valueOf(delta) + " "
								+ String.valueOf(TPN));

						maps[i] = map;
					}
				}

				continue;
			}

			if (i != no && delta * 100.0 / range > TP && TP > 0.0) {
				if (info)
					System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(delta));
					map.put("mask_reason", "Current point's delta distance > Pref.TP:"
							+ String.valueOf(delta * 100.0 / range) + " " + String.valueOf(TP));

					maps[i] = map;
				}

				continue;
			} else if (vFlag && t < THETA && delta * 100.0 / range > TPV && THETA > 0.0 && TPV > 0.0) {
				if (info)
					System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"
							+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(t));
					map.put("mask_reason", "Current point's delta angle > Pref.TPV:" + String.valueOf(t) + " "
							+ String.valueOf(THETA));

					maps[i] = map;
				}
			}

			if (ds != null) {
				if (delta > 499)
					delta = 499;
				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;

		}

		if (HillConstants.bellMask) {
			// check bell curve
			boolean bellFlag = checkBellShape(xs, ys, sd);

			if (bellFlag) {
				maxY = 0.0;
				for (int i = ys.length / 2; i < ys.length; i++) {
					if (Math.abs(ys[i]) > maxY)
						maxY = Math.abs(ys[i]);
				}

				// assign masked point
				for (int j = 0; j < list.size(); j++) {
					String line = (String) list.get(j);
					int k = line.indexOf(":");

					int l = Integer.parseInt(line.substring(0, k));
					flags[l] = false;
				}

				for (int i = ys.length - 1; i >= ys.length / 2; i--) {
					// unmask if masked
					if (Math.abs(Math.abs(ys[i]) - maxY) < 1.0e-6) {
						flags[i] = true;
						maps[i] = null;

						break;
					}

					if (flags[i] == false)
						continue;

					if (Math.abs(ys[i]) < maxY - sd) // *0.95)
					{
						list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
						flags[i] = false;

						if (maps != null) {
							HashMap map = new HashMap();

							map.put("sample_conc", String.valueOf(xs[i]));
							map.put("sample_resp", String.valueOf(ys[i]));
							map.put("mask_value", String.valueOf(ys[i]));
							map.put("mask_reason", "y<maxY:" + String.valueOf(ys[i]) + " " + String.valueOf(PMAX));

							maps[i] = map;
						}
					}
				}

				return flags;
			}
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++) {
			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;

			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++) {
				String line = (String) list.get(j);
				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd) {
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);

		}

		for (int i = 0; i < xs.length; i++) {
			if (flags[i] && maps != null)
				maps[i] = null;
		}

		return flags;
	}

	public static boolean checkBellShape(double[] xs, double[] ys1, double sd) {
		double[] ys = new double[ys1.length];
		for (int i = 0; i < ys1.length; i++)
			ys[i] = Math.abs(ys1[i]);

		// check bell curve
		double maxY = 0.0;
		int pos = 0;
		for (int i = ys.length / 2; i < ys.length; i++) {
			if (ys[i] > maxY) {
				maxY = ys[i];

				pos = i;
			}
		}

		// first check a supporting point
		if (pos > 1 && ys[pos] > ys[pos - 1] && ys[pos - 1] > ys[pos - 2] && ys[pos - 1] > sd)
			return true;

		// second check the data points after
		SimpleRegression sr = new SimpleRegression();
		int no = 0;
		for (int j = pos + 1; j < ys.length; j++) {
			sr.addData(xs[j], ys[j]);
			no++;
		}

		double r = sr.getR();
		double s = sr.getSlope();
		double c = sr.getIntercept();
		double p = 1.0;

		try {
			if (no > 2)
				p = sr.getSignificance();
		} catch (Exception ex) {
			ex.printStackTrace();
		}

		if (p < 0.05 && s > 0.0)
			return false;

		return true;
	}

	public static boolean[] maskDif3(double[] xs, double[] ys, int[][] ds, double[] ms, HashMap[] maps) {
		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		// calc mean and std
		int no = 0;

		double mean = 0.0;
		for (int i = 0; i < xs.length; i++) {
			mean += ys[i];
			no++;
		}
		mean /= no;

		double dev = 0.0;
		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - mean;
			dev += delta * delta;
		}
		dev = Math.sqrt(dev / no);

		// check first point

		// if (false && P0 >0.0 && Math.abs(ys[0]) > P0)
		if (PMAX > 0.0 && Math.abs(ys[0]) > PMAX) {
			if (info)
				System.out.println("Abs(y0)>Pref.PMAX");
			// System.out.println("Abs(y0)>Pref.P0");

			flags[0] = false;

			if (maps != null) {
				HashMap map = new HashMap();

				map.put("sample_conc", String.valueOf(xs[0]));
				map.put("sample_resp", String.valueOf(ys[0]));
				map.put("mask_value", String.valueOf(ys[0]));
				map.put("mask_reason", "abs(y0)>Pref.PMAX:" + String.valueOf(ys[0]) + " " + String.valueOf(PMAX));
				maps[0] = map;
			}
		}

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++) {
			// check against absolute value
			if (PMAX > 0.0 && Math.abs(ys[i]) > PMAX) {
				if (info)
					System.out.println("Abs(y)>Pref.PMAX");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(ys[i]));
					map.put("mask_reason", "abs(y)>Pref.PMAX:" + String.valueOf(ys[i]) + " " + String.valueOf(PMAX));

					maps[i] = map;
				}

				continue;
			}

			// check against deviation from mean
			double ratio = Math.abs(ys[i] - mean) / dev;
			if (TR > 0.0 && ratio > TR) {
				if (info)
					System.out.println("(resp-mean)/std > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(ratio));
					map.put("mask_reason", "(resp-mean)/std > Pref.TR:" + String.valueOf(ratio) + " "
							+ String.valueOf(TR));

					maps[i] = map;
				}

				continue;
			}

			// shape like v
			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
				vFlag = true;

			// check against angle and delta-the difference from the
			// extrapolated line
			double t1 = 0.0;
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;
			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;
			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta1 = 10000000.0;
			double delta2 = 10000000.0;

			if (i == no) {
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2])
						/ (xs[no - 1] - xs[no - 2]) - ys[no - 2]);
			} else {
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
						- ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;

			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			if (debug)
				System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no) {
				if (delta > TPN && TPN > 0.0) {
					if (info)
						System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

					if (maps != null) {
						HashMap map = new HashMap();

						map.put("sample_conc", String.valueOf(xs[i]));
						map.put("sample_resp", String.valueOf(ys[i]));
						map.put("mask_value", String.valueOf(delta));
						map.put("mask_reason", "Last point's delta distance > Pref.TPN:" + String.valueOf(delta) + " "
								+ String.valueOf(TPN));
						maps[i] = map;
					}
				}

				continue;
			}

			if (i != no && delta > TP && TP > 0.0) {
				if (info)
					System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(delta));
					map.put("mask_reason", "Current point's delta distance > Pref.TP:" + String.valueOf(delta) + " "
							+ String.valueOf(TP));

					maps[i] = map;
				}

				continue;
			} else if (vFlag && t < THETA && delta > TPV && THETA > 0.0 && TPV > 0.0) {
				if (info)
					System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"
							+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) {
					HashMap map = new HashMap();

					map.put("sample_conc", String.valueOf(xs[i]));
					map.put("sample_resp", String.valueOf(ys[i]));
					map.put("mask_value", String.valueOf(t));
					map.put("mask_reason", "Current point's delta angle > Pref.TPV:" + String.valueOf(t) + " "
							+ String.valueOf(THETA));

					maps[i] = map;
				}
			}

			if (ds != null) {
				if (delta > 499)
					delta = 499;

				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++) {
			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;
			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++) {
				String line = (String) list.get(j);
				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd) {
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}

		for (int i = 0; i < xs.length; i++) {
			if (flags[i] && maps != null)
				maps[i] = null;
		}

		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);

		return flags;
	}

	public static boolean[] maskDif1(double[] xs, double[] ys, int[][] ds, double[] ms) {
		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		int no = 0;

		double mean = 0.0;
		for (int i = 0; i < xs.length; i++) {
			mean += ys[i];
			no++;
		}
		mean /= no;

		double dev = 0.0;
		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - mean;
			dev += delta * delta;
		}
		dev = Math.sqrt(dev / no);

		if (P0 > 0.0 && Math.abs(ys[0]) > P0) {
			if (info) {
				System.out.println("ABS of y0 >Pref.P0");
			}

			flags[0] = false;
		}

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		for (int i = 1; i < xs.length; i++) {
			if (Math.abs(ys[i]) > PMAX) {
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				continue;
			}

			// based upon deviation from mean
			double ratio = Math.abs(ys[i] - mean) / dev;
			if (debug)
				System.out.println(i + ":" + ratio);
			if (ratio > TR) {
				if (info)
					System.out.println("Deviation from mean > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				continue;
			}

			// based upon delta
			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta2 = 10000000.0;
			double delta1 = 10000000.0;

			// shape like v
			boolean vFlag = false;
			double t1 = 0.0;
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;

			if (i < no) {
				if ((ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
					vFlag = true;

				t2 = calcTheta(xs, ys, i - 1, i, i + 1);
			}

			double t = 0.0;
			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			if (i == no) {
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2])
						/ (xs[no - 1] - xs[no - 2]) - ys[no - 2]);
			} else {
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1])
						- ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;
			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			if (debug)
				System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (i == no || i == no - 1) {
				if (ys[i] - ys[i - 1] > 0.0 && delta > TPV) {
					if (info)
						System.out.println(delta + " First and last point's delta > Pref.TPV");
					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				} else if (delta > TPN) {
					if (info)
						System.out.println("Current point's delta > Pref.TPN");
					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				}
			} else if (delta > TP) {
				if (info)
					System.out.println("Current point's delta > Pref.TP");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
			} else if (t < THETA && delta > TPV) {
				if (info)
					System.out.println("Current point's delta > Pref.TPV2");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
			}

			if (ds != null) {
				if (delta > 499)
					delta = 499;
				ds[i][(int) delta]++;
			}
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++) {
			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;
			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++) {
				String line = (String) list.get(j);
				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd) {
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}

		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);

		return flags;
	}

	public static boolean[] maskR2(double[] xs, double[] ys, int[][] ds, double[] ms) {
		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		// always mask bad first point
		if (P0 > 0.0 && Math.abs(ys[0]) > P0)
			flags[0] = false;

		// if r2 is already good, do not mask any more
		double[] values = hillFit(PI_MAX, PS_MIN, flags, null, xs, ys);

		if (values == null || values[3] > R2)
			return flags;

		// how many points allowed for masking
		int no = xs.length / MASK;

		if (debug)
			System.out.println(no);

		boolean breakFlag = false;

		for (int j = 0; j < no; j++) {
			boolean[] flags2 = new boolean[xs.length];

			for (int i = 0; i < xs.length; i++)
				flags2[i] = flags[i];

			int maxI = -1;
			double maxr2 = 0.0;

			for (int i = 1; i < xs.length; i++) {
				if (flags[i] == false)
					continue;

				flags2[i] = false;

				values = hillFit(PI_MAX, PS_MIN, flags2, null, xs, ys);

				if (values != null)
					System.out.println(i + ":" + values[3]);

				if (values == null) {
					maxI = i;
					breakFlag = true;
					break;
				} else if (values[3] > maxr2) {
					maxr2 = values[3];
					maxI = i;
				}

			}

			flags[maxI] = false;

			if (breakFlag)
				break;
		}

		return flags;
	}

	public static double dist(double xk, double yk, double xl, double yl) {
		double dy = 0.1 * (yk - yl);
		double dx = xk - xl;

		return Math.sqrt(dx * dx + dy * dy);
	}

	public static double dist(double[] xs, double[] ys, int k, int l) {
		double dy = 0.1 * (ys[k] - ys[l]);
		double dx = xs[k] - xs[l];

		return Math.sqrt(dx * dx + dy * dy);
	}

	public static double calcTheta(double[] xs, double[] ys, int k, int l, int m) {
		double dkl = dist(xs, ys, k, l);
		double dlm = dist(xs, ys, l, m);
		double dkm = dist(xs, ys, k, m);

		double cos = (dkl * dkl + dlm * dlm - dkm * dkm) / (2 * dkl * dlm);

		return 57.0 * Math.acos(cos);

	}

	public static double calcTheta(double xk, double yk, double xl, double yl, double xm, double ym) {
		double dkl = dist(xk, yk, xl, yl);
		double dlm = dist(xl, yl, xm, ym);
		double dkm = dist(xk, yk, xm, ym);

		double cos = (dkl * dkl + dlm * dlm - dkm * dkm) / (2 * dkl * dlm);

		return 57.0 * Math.acos(cos);
	}

	public static boolean[] maskAngle(double[] xs, double[] ys, int[][] ds, double[] ms) {
		int n = xs.length - 1;

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (P0 > 0.0 && Math.abs(ys[0]) > P0)
			flags[0] = false;

		for (int i = 1; i < n; i++) {
			double theta = 0.0;

			if (i == 1)
				theta = calcTheta(xs[i - 1], 0.0, xs[i], ys[i], xs[i + 1], ys[i + 1]);
			else
				theta = calcTheta(xs, ys, i - 1, i, i + 1);

			ds[i][(int) theta]++;

		}

		return flags;

	}

	public static boolean[] maskLine(double[] xs, double[] ys) {
		boolean[] flag1 = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flag1[i] = true;

		double[] fitValues = HtsUtil.lineFit(xs, ys);

		if (debug)
			System.out.println(fitValues[0] + ":" + fitValues[1] + ":" + fitValues[2]);

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - fitValues[1] - fitValues[0] * xs[i];
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / xs.length);

		boolean[] flag2 = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++) {
			double y0 = fitValues[1] - fitValues[0] * xs[i];

			if (Math.abs(ys[i] - y0) > 3 * dev)
				flag2[i] = false;
			else
				flag2[i] = true;

			double ratio = Math.abs(ys[i] - y0) / dev;

			if (debug)
				System.out.println(i + ":" + y0 + ":" + dev + ":" + (ys[i] - y0) + ":" + ratio + ":" + flag2[i]);

		}

		return flag2;

	}

	public static boolean checkRange(boolean[] flags, double[] ys) {
		double y_min = 1000000.0;

		double y_max = -1000000.0;

		for (int i = 0; i < ys.length; i++)

			if (flags == null || flags[i]) {
				if (ys[i] < y_min)
					y_min = ys[i];

				if (ys[i] > y_max)
					y_max = ys[i];
			}

		if (y_max - y_min < MIN_Y_RANGE)
			return false;

		return true;

	}

	public static boolean checkRange(double ymin, double ymax, boolean[] flags, double[] ys) {
		double y_min = Double.MAX_VALUE;

		double y_max = -Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)

			if (flags == null || flags[i])

			{

				if (ys[i] < y_min)

					y_min = ys[i];

				if (ys[i] > y_max)

					y_max = ys[i];

			}

		double range = y_max - y_min;

		if (y_min > HillConstants.SUPER_Y || y_max < -HillConstants.SUPER_Y) {
			return true;
		}

		// if ( 100.0*range/Math.max(Math.abs(y_min), Math.abs(y_max)) <
		// MIN_Y_RANGE)
		if (range < MIN_Y_RANGE) {
			// System.out.println("range check: "+(100.0*range/Math.max(Math.abs(y_min),
			// Math.abs(y_max))));
			return false;
		}

		return true;

	}

	public static boolean isOnePoint(double[] ys) {
		int no = 0;

		for (int i = 0; i < ys.length; i++) {
			if (Math.abs(ys[i]) < HillConstants.Y0)
				no++;
		}

		if (no < ys.length - 1)
			return false;

		if (ys.length < 3)
			return false;

		if (Math.abs(ys[ys.length - 1]) < HillConstants.Y0 && Math.abs(ys[ys.length - 2]) < HillConstants.Y0)
			return true;

		return false;
	}

	public static ArrayList getMaskList(boolean[] flags, double[] xs, double[] ys) {
		ArrayList list = new ArrayList();

		for (int i = 0; i < flags.length; i++) {
			if (flags[i] == false)
				list.add("0-" + String.valueOf(i));
		}

		return list;
	}

	public static double[] lineFit(boolean flags[], double[] xs, double ys[]) {
		int no = xs.length;

		double xm = 0;
		double ym = 0;

		int nn = 0;

		for (int i = 0; i < no; i++)

			if (flags == null || flags[i]) {
				xm += xs[i];
				ym += ys[i];
				nn++;
			}

		xm /= nn;
		ym /= nn;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;

		for (int i = 0; i < no; i++)
			if (flags == null || flags[i]) {
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

	public static double[] partialHillFit(String[] labels, String[] values) throws Exception {
		ArrayList listx = new ArrayList();
		ArrayList listy = new ArrayList();

		double maxY = -1000000.0;
		int maxI = 0;

		int nn = 0;

		for (int i = 0; i < labels.length; i++) {

			if (labels[i] == null || labels[i].equals("null"))

				throw new Exception("Missing labels");

			if (values[i] != null && values[i].length() > 2) {
				double conc = HtsUtil.parseConcentrationMolar(labels[i]);

				double ty = Double.parseDouble(values[i]);

				if (ty > maxY) {
					maxI = i;
					maxY = ty;
				}

				listx.add(String.valueOf(Math.log10(conc)));
				;

				listy.add(values[i]);

			} else
				nn++;

		}

		int no = maxI - nn + 1;

		// require minimum 4 data points
		if (no < 4)
			return null;

		double[] x = new double[no];
		double[] y = new double[no];

		for (int i = 0; i < no; i++) {
			x[i] = Double.parseDouble((String) listx.get(i));
			y[i] = Double.parseDouble((String) listy.get(i));
		}

		double[] fitValues = partialHillFit(maxY, x, y);

		return fitValues;
	}

	// 2p fit
	public static double[] partialHillFit(double maxY, double[] xs, double[] ys) {
		if (checkRange(null, ys) == false)
			return null;

		double yinf_min = 0.0;
		double slope_min = 1.0;
		double x05_min = 0.0;
		double dev_min = 1000000.0;

		double yl = 0.0;
		double yr = HillConstants.Y0_INF_COEF * maxY;
		double yd = 1.0;

		for (double yinf = yl; yinf < yr; yinf += yd) {
			for (double x05 = -10.0; x05 < -4.0; x05 += 0.1) {
				double dev = calcHillDeviation(x05, 0.0, yinf, 1.0, null, null, xs, ys);

				// System.out.println(yinf+":"+x05+":"+slope+":"+dev);

				if (dev < dev_min) {
					dev_min = dev;
					x05_min = x05;
					yinf_min = yinf;
				}
			}
		}

		double dev_const = calcConstantDeviation(null, null, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		// double r2 = 1.0-Math.sqrt(dev_min/dev_const);

		// System.out.println(r2+":"+R2);

		if (r2 < R2)

			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, 0.0 };

		return values;

	}

	public static double[] hillFit(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		if (HillConstants.bellMask)
			return hillFit(ymin, ymax, flags, ws, xs, ys, true, true);
		else
			return hillFit(ymin, ymax, flags, ws, xs, ys, false, true);
	}

	public static double[] hillFit(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys,
			boolean flag2, boolean p4Fit) {
		return hillFit(ymin, ymax, flags, ws, xs, ys, flag2, String.valueOf(p4Fit));

	}

	public static double[] hillFit(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys,
			boolean flag2, String option) {
		double[] values = null;

		if (option == null) {
			if (HillConstants.P4_FIT)
				option = "true";
			else
				option = "false";
		}

		if (option.equals("true") || option.equals("false"))
			values = hillFitFast(ymin, ymax, flags, ws, xs, ys, Boolean.valueOf(option));
		else
			values = hillFitFast(ymin, ymax, flags, ws, xs, ys, option);

		if (values == null)
			return null;

		if (flag2 || HillConstants.bellMask)
			return values;

		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;

		for (int i = 0; i < xs.length; i++) {
			max = Math.max(max, ys[i]);
			min = Math.min(min, ys[i]);
		}
		double range = max - min;

		double y0 = values[6];
		double x05 = values[0];
		double yinf = values[1];
		double slope = values[2];

		double[] yfit = calcHillFitCurve(x05, y0, yinf, slope, xs);

		boolean flag = false;

		// first unmask good points
		int maskNo = 0;

		for (int i = 1; i < ys.length; i++) {
			if (flags[i] == false && Math.abs(ys[i] - yfit[i]) * 100.0 / range < TPV) {
				flag = true;
				flags[i] = true;
			}

			if (flags[i] == false)
				maskNo++;
		}

		ArrayList list = new ArrayList();

		if (maskNo < xs.length / MASK) {
			for (int i = 0; i < ys.length - 1; i++) {
				if (i > 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TPV) {
					// flag = true;

					list.add(String.valueOf(i));

					// flags[i]= false;

				} else if (i == 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TP0) {
					// flag = true;

					list.add(String.valueOf(i));

					// flags[i]= false;
				}
			}
		}

		if (list.size() > 0) {
			flag = true;

			for (int i = maskNo; i < xs.length / MASK && (i - maskNo) < list.size(); i++) {
				int k = Integer.parseInt((String) list.get(i - maskNo));

				flags[k] = false;
			}
		}

		if (flag) {
			if (option.equals("true") || option.equals("false"))
				return hillFitFast(ymin, ymax, flags, ws, xs, ys, Boolean.valueOf(option));
			else
				return hillFitFast(ymin, ymax, flags, ws, xs, ys, option);
		}

		return values;
	}

	// biphase fit
	public static double[] hillFitFast(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs,
			double[] ys, String option) {
		// mono phase fit
		if (option.toLowerCase().startsWith("mono"))
			return hillFitFast(ymin, ymax, flags, ws, xs, ys, true);

		if (checkRange(ymin, ymax, flags, ys) == false)
			return null;

		// biphase fit
		double xl = -10.0;
		double xr = -2.0;
		double xd = 0.5;

		double y0l = -150.0;
		double y0r = 150.0;
		double y0d = 5.0;

		double yil = -150.0;
		double yir = 150.0;
		double yid = 5.0;

		double s1l = MIN_SLOPE;
		double s1r = MAX_SLOPE;
		double s1d = 1.0;

		double fl = 0.0;
		double fr = 1.0;
		double fd = 0.1;

		double s2l = MIN_SLOPE;
		double s2r = MAX_SLOPE;
		double s2d = 0.5;

		// min and max
		double max_y = -Double.MAX_VALUE;
		double min_y = Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)
			if (flags == null || flags[i]) {
				max_y = Math.max(max_y, ys[i]);
				min_y = Math.min(min_y, ys[i]);
			}

		double[] fitValues = lineFit(flags, xs, ys);

		if (fitValues[0] > 0) {
			if (min_y < 0) {
				y0l = HillConstants.Y0_INF_COEF * min_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = HillConstants.Y0_INF_COEF * min_y;
			}

			if (max_y > 0.0) {
				yil = 0.0;
				yir = HillConstants.Y0_INF_COEF * max_y;
			} else {
				yir = 0.0;
				yil = HillConstants.Y0_INF_COEF * max_y;
			}
		} else {
			if (max_y < 0) {
				y0l = HillConstants.Y0_INF_COEF * max_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = HillConstants.Y0_INF_COEF * max_y;
			}

			if (min_y < 0.0) {
				yil = HillConstants.Y0_INF_COEF * min_y;
				yir = 0;
			} else {
				yir = HillConstants.Y0_INF_COEF * min_y;
				yil = 0;
			}
		}

		double delta_y = 0.05 * (max_y - min_y);
		if (delta_y < 5.0)
			delta_y = 5.0;

		y0d = delta_y;
		yid = delta_y;

		// estimate
		fitValues = hillFit2(y0l, y0r, y0d, yil, yir, yid, fl, fr, fd, xl, xr, xd, flags, ws, xs, ys);

		if (fitValues == null)
			return null;

		// double[] values =
		// {y0_min,yi_min,f_min,e1_min,s1_min,e2_min,s2_min,r2,dev_const,dev_min};
		fl = fitValues[2] - 0.1;
		fr = fitValues[2] + 0.1;
		fd = 0.05;

		double x1l = fitValues[3] - 1.0;
		double x1r = fitValues[3] + 1.0;
		double x1d = HillConstants.EC50_DELTA;

		double x2l = fitValues[5] - 1.0;
		double x2r = fitValues[5] + 1.0;
		double x2d = HillConstants.EC50_DELTA;

		double y0 = fitValues[0];
		double yi = fitValues[1];

		s1l = 0.2;
		s1r = 4.0;
		s1d = 0.2;

		s2l = 0.2;
		s2r = 4.0;
		s2d = 0.2;

		// refine slope and ec50
		fitValues = hillFit2(y0, yi, fl, fr, fd, x1l, x1r, x1d, x2l, x2r, x2d, s1l, s1r, s1d, s2l, s2r, s2d, flags, ws,
				xs, ys);

		// refine y0 and yi
		y0l = fitValues[0] - 5.0;
		y0r = fitValues[0] + 5.0;
		y0d = 1.0;

		yil = fitValues[1] - 5.0;
		yir = fitValues[1] + 5.0;
		yid = 1.0;

		fl = fitValues[2] - 0.05;
		fr = fitValues[2] + 0.05;
		fd = 0.01;

		double x1 = fitValues[3];
		double x2 = fitValues[5];
		double s1 = fitValues[4];
		double s2 = fitValues[6];

		fitValues = hillFit2(y0l, y0r, y0d, yil, yir, yid, fl, fr, fd, x1, x2, s1, s2, flags, ws, xs, ys);

		return fitValues;
	}

	public static double[] hillFitFast(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs,
			double[] ys, boolean p4Fit) {
		if (checkRange(ymin, ymax, flags, ys) == false)
			return null;

		double xl = -10.0;
		double xr = -2.0;
		double xd = 0.5;

		double y0l = -150.0;
		double y0r = 150.0;
		double y0d = 5.0;

		double yl = -150.0;
		double yr = 150.0;
		double yd = 5.0;

		double sl = MIN_SLOPE;
		double sr = MAX_SLOPE;
		double sd = 0.2;

		double[] fitValues = lineFit(flags, xs, ys);

		double max_y = -Double.MAX_VALUE;
		double min_y = Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)
			if (flags[i]) {
				if (ys[i] > max_y)
					max_y = ys[i];

				if (ys[i] < min_y)
					min_y = ys[i];
			}

		double delta_y = 0.05 * (max_y - min_y);
		if (delta_y < 2.0)
			delta_y = 2.0;

		if (fitValues[0] > 0) {
			if (min_y < 0) {
				y0l = HillConstants.Y0_INF_COEF * min_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = HillConstants.Y0_INF_COEF * min_y;
			}

			if (max_y > 0.0) {
				yl = 0.0;
				yr = HillConstants.Y0_INF_COEF * max_y;
			} else {
				yr = 0.0;
				yl = HillConstants.Y0_INF_COEF * max_y;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		} else {
			if (max_y < 0) {
				y0l = HillConstants.Y0_INF_COEF * max_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = HillConstants.Y0_INF_COEF * max_y;
			}

			if (min_y < 0.0) {
				yl = HillConstants.Y0_INF_COEF * min_y;
				yr = 0;
			} else {
				yr = HillConstants.Y0_INF_COEF * min_y;
				yl = 0;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		}

		y0d = delta_y;
		yd = delta_y;

		if (info)
			System.out.println(y0l + ":" + y0l + ":" + yl + ":" + yr + ":" + xl + ":" + xr);

		// System.out.println(p4Fit+":"+min_y+":"+max_y+":"+delta_y);

		if (p4Fit)
			fitValues = hillFit(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys);
		else
			fitValues = hillFitY0(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, flags, ws, xs, ys);

		if (fitValues == null)
			return null;

		xl = fitValues[0] - 1.0;
		xr = fitValues[0] + 1.0;

		xd = HillConstants.EC50_DELTA;

		y0l = fitValues[6] - 2 * delta_y;
		y0r = fitValues[6] + 2 * delta_y;

		yl = fitValues[1] - 2 * delta_y;
		yr = fitValues[1] + 2 * delta_y;

		delta_y = 0.1 * delta_y;
		if (delta_y < 0.5)
			delta_y = 0.5;

		y0d = delta_y;
		yd = delta_y;

		if (fitValues[2] - 0.5 < 0.1)
			sl = 0.1;
		else
			sl = fitValues[2] - 0.5;

		sr = fitValues[2] + 0.5;
		sd = 0.1;

		if (p4Fit)
			fitValues = hillFit(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys);
		else
			fitValues = hillFitY0(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, flags, ws, xs, ys);

		return fitValues;
	}

	// 2p fit
	public static double[] hillFit2P(double maxY, double[] xs, double[] ys) {
		if (checkRange(null, ys) == false)
			return null;

		double yinf_min = 0.0;
		double slope_min = 1.0;
		double x05_min = 0.0;
		double dev_min = 1000000.0;

		double yl = 0.0;
		double yr = HillConstants.Y0_INF_COEF * maxY;
		double yd = 1.0;

		for (double yinf = yl; yinf < yr; yinf += yd) {
			for (double x05 = -10.0; x05 < -4.0; x05 += 0.1) {
				double dev = calcHillDeviation(x05, 0.0, yinf, 1.0, null, null, xs, ys);

				if (dev < dev_min) {
					dev_min = dev;
					x05_min = x05;
					yinf_min = yinf;
				}
			}
		}

		double dev_const = calcConstantDeviation(null, null, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, 0.0 };

		return values;

	}

	// 3p fit with y0 fixed
	public static double[] hillFitSlope(double yl, double yr, double yd, double xl, double xr, double xd, double sl,
			double sr, double sd, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 1000000.0;

		double delta = sd;

		int no = 0;

		for (double yinf = yl; yinf <= yr; yinf += yd) {
			for (double x05 = xl; x05 <= xr; x05 += xd) {
				for (double slope = sl; slope <= sr; slope += delta) {
					if (Math.abs(slope) < MIN_SLOPE)
						continue;

					double dev = calcHillDeviation(x05, 0.0, yinf, slope, flags, ws, xs, ys);

					if (dev < dev_min) {
						dev_min = dev;
						slope_min = slope;
						x05_min = x05;
						yinf_min = yinf;
					}

					delta = Math.abs(slope) * 0.1;

					if (delta < 0.1)
						delta = 0.1;

					no++;
				}
			}

		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, 0.0 };

		return values;

	}

	// 3p fit with slope fixed
	public static double[] hillFitY0(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 100000000.0;

		double slope = 1.0;

		int no = 0;

		for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
			for (double yinf = yl; yinf <= yr; yinf += yd) {
				for (double x05 = xl; x05 <= xr; x05 += xd) {
					double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);

					if (dev < dev_min) {
						dev_min = dev;
						slope_min = slope;
						x05_min = x05;
						yinf_min = yinf;
						y0_min = y0;
					}

					no++;
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;
		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope, r2, dev_const, dev_min, y0_min };

		return values;
	}

	public static double calcSlope(double x05, double y0, double yinf, boolean[] flags, double[] xs, double[] ys) {
		double[] x = new double[xs.length];
		double[] y = new double[xs.length];

		int no = 0;

		ArrayList list = new ArrayList();

		for (int i = 0; i < xs.length; i++) {
			if (flags != null && flags[i] == false)
				continue;

			if (Math.abs(ys[i] - y0) < 0.01)
				continue;

			double td1 = (yinf - ys[i]) / (ys[i] - y0);
			if (td1 < 0.0)
				continue;

			td1 = Math.log10(td1);

			double td2 = x05 - xs[i];

			if (no > 0 && Math.abs(td1 - y[no - 1]) < 0.01) {
				if (no > 1) {
					double slope = calcSlope(no, x, y);
					list.add(String.valueOf(slope));
				}

				no = 0;
			}

			// System.out.printf("%5d %8.3f %8.3f\n",i,td1,td2);

			x[no] = td2;
			y[no] = td1;

			no++;
		}

		if (no > 1) {
			double slope = calcSlope(no, x, y);
			list.add(String.valueOf(slope));
		}

		if (list.size() == 0)
			return -1;

		return calcMedian(list);
	}

	public static double calcSlope(int no, double[] x, double[] y) {
		double xm = 0;
		double ym = 0;

		int k = 0;
		for (int i = 0; i < no; i++) {
			xm += x[i];
			ym += y[i];

			k++;
		}

		xm /= k;
		ym /= k;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;
		for (int i = 0; i < no; i++) {
			lxx += (x[i] - xm) * (x[i] - xm);
			// lyy += (y[i]-ym)*(y[i]-ym);
			lxy += (x[i] - xm) * (y[i] - ym);
		}

		double b = lxy / lxx;
		// double a = ym-b*xm;
		// double r = lxy/Math.sqrt(lxx*lyy);

		return b;
	}

	public static double calcMedian(ArrayList list) {
		double[] values = new double[list.size()];

		for (int k = 0; k < list.size(); k++) {
			values[k] = Double.parseDouble((String) list.get(k));
		}

		Arrays.sort(values);

		return values[values.length / 2];
	}

	// 4p fit
	public static double[] hillFitFast(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, double sl, double sr, double sd, boolean[] flags, double[] ws, double[] xs,
			double[] ys) {
		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 100000000.0;
		double delta = sd;

		int no = 0;

		for (double x05 = xl; x05 <= xr; x05 += xd) {
			for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
				for (double yinf = yl; yinf <= yr; yinf += yd) {
					double s0 = calcSlope(x05, y0, yinf, flags, xs, ys);

					if (s0 < MIN_SLOPE)
						s0 = MIN_SLOPE;

					if (s0 > MAX_SLOPE)
						s0 = MAX_SLOPE;

					// at s0
					double dev0 = calcHillDeviation(x05, y0, yinf, s0, flags, ws, xs, ys);
					if (dev0 < dev_min) {
						dev_min = dev0;
						slope_min = s0;
						x05_min = x05;
						yinf_min = yinf;
						y0_min = y0;
					}

				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;
		// if (r2 < R2)
		// return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };

		return values;
	}

	public static void setFastFlag(boolean flag) {
		fastFlag = flag;
	}

	// 4p fit
	public static double[] hillFit(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, double sl, double sr, double sd, boolean[] flags, double[] ws, double[] xs,
			double[] ys) {
		if (fastFlag)
			return hillFitFast(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys);

		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 100000000.0;

		double delta = sd;

		int no = 0;

		// any parameter fixed?
		if (HillConstants.FIXED_Y0 != null && HillConstants.FIXED_Y0.length() > 0) {
			double td = Double.parseDouble(HillConstants.FIXED_Y0);
			y0l = td;
			y0r = td;
		}

		if (HillConstants.FIXED_YINF != null && HillConstants.FIXED_YINF.length() > 0) {
			double td = Double.parseDouble(HillConstants.FIXED_YINF);
			yl = td;
			yr = td;
		}

		if (HillConstants.FIXED_SLOPE != null && HillConstants.FIXED_SLOPE.length() > 0) {
			double td = Double.parseDouble(HillConstants.FIXED_SLOPE);
			sl = td;
			sr = td;
		}

		for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
			for (double yinf = yl; yinf <= yr; yinf += yd) {
				for (double x05 = xl; x05 <= xr; x05 += xd) {
					for (double slope = sl; slope <= sr; slope += delta) {
						if (Math.abs(slope) < MIN_SLOPE)
							continue;

						double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);

						if (dev < dev_min) {
							dev_min = dev;
							slope_min = slope;
							x05_min = x05;
							yinf_min = yinf;
							y0_min = y0;
						}

						delta = Math.abs(slope) * 0.1;
						if (delta < 0.1)
							delta = 0.1;

						no++;
					}
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;
		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };

		return values;
	}

	// biphase fit
	public static double[] hillFit2(double y0l, double y0r, double y0d, double yil, double yir, double yid, double fl,
			double fr, double fd, double xl, double xr, double xd, boolean[] flags, double[] ws, double[] xs,
			double[] ys) {
		double y0_min = 0.0;
		double yi_min = 0.0;
		double f_min = 0.0;
		double e1_min = 0.0;
		double e2_min = 0.0;
		double dev_min = Double.MAX_VALUE;

		int no = 0;

		double s1 = 1.0;
		double s2 = 1.0;

		for (double y0 = y0l; y0 <= y0r; y0 += y0d)
			for (double yi = yil; yi <= yir; yi += yid)
				for (double f = fl; f <= fr; f += fd)
					for (double x1 = xl; x1 <= xr; x1 += xd)
						for (double x2 = x1 + xd; x2 <= xr; x2 += xd) {
							double dev = calcHillDeviation(y0, yi, f, x1, s1, x2, s2, flags, ws, xs, ys);

							if (dev < dev_min) {
								dev_min = dev;
								y0_min = y0;
								yi_min = yi;

								e1_min = x1;
								e2_min = x2;
								f_min = f;
							}

							no++;
						}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		double[] values = { y0_min, yi_min, f_min, e1_min, 1.0, e2_min, 1.0, r2, dev_const, dev_min };

		if (r2 < R2)
			return null;

		return values;
	}

	// biphase fit
	public static double[] hillFit2(double y0l, double y0r, double y0d, double yil, double yir, double yid, double fl,
			double fr, double fd, double x1, double x2, double s1, double s2, boolean[] flags, double[] ws,
			double[] xs, double[] ys) {
		double y0_min = 0.0;
		double yi_min = 0.0;
		double f_min = 0.0;
		double s1_min = 0.0;
		double s2_min = 0.0;
		double e1_min = 0.0;
		double e2_min = 0.0;
		double dev_min = Double.MAX_VALUE;

		int no = 0;

		for (double y0 = y0l; y0 <= y0r; y0 += y0d)
			for (double yi = yil; yi <= yir; yi += yid)
				for (double f = fl; f <= fr; f += fd) {
					double dev = calcHillDeviation(y0, yi, f, x1, s1, x2, s2, flags, ws, xs, ys);

					if (dev < dev_min) {
						dev_min = dev;
						y0_min = y0;
						yi_min = yi;
						s1_min = s1;
						s2_min = s2;
						e1_min = x1;
						e2_min = x2;
						f_min = f;
					}

					no++;
				}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		double[] values = { y0_min, yi_min, f_min, e1_min, s1_min, e2_min, s2_min, r2, dev_const, dev_min };

		if (r2 < R2)
			return null;

		return values;
	}

	// biphase fit
	public static double[] hillFit2(double y0, double yi, double fl, double fr, double fd, double x1l, double x1r,
			double x1d, double x2l, double x2r, double x2d, double s1l, double s1r, double s1d, double s2l, double s2r,
			double s2d, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double y0_min = 0.0;
		double yi_min = 0.0;
		double f_min = 0.0;
		double s1_min = 0.0;
		double s2_min = 0.0;
		double e1_min = 0.0;
		double e2_min = 0.0;
		double dev_min = Double.MAX_VALUE;

		double delta1 = s1d;
		double delta2 = s2d;

		int no = 0;

		// for (double y0 = y0l; y0<=y0r; y0+=y0d)
		// for (double yi = yil; yi<=yir; yi+=yid)
		for (double f = fl; f <= fr; f += fd)
			for (double s1 = s1l; s1 <= s1r; s1 += delta1)
				for (double s2 = s2l; s2 <= s2r; s2 += delta2)
					for (double x1 = x1l; x1 <= x1r; x1 += x1d)
						for (double x2 = x2l; x2 <= x2r; x2 += x2d) {
							double dev = calcHillDeviation(y0, yi, f, x1, s1, x2, s2, flags, ws, xs, ys);

							if (dev < dev_min) {
								dev_min = dev;
								y0_min = y0;
								yi_min = yi;
								s1_min = s1;
								s2_min = s2;
								e1_min = x1;
								e2_min = x2;
								f_min = f;
							}

							delta1 = s1 * 0.1;
							if (delta1 < 0.2)
								delta1 = 0.2;

							delta2 = s2 * 0.1;
							if (delta2 < 0.2)
								delta2 = 0.2;

							// if (no % 10000 == 0)
							// System.out.println(no+"\t"+y0+"\t"+yi+"\t"+f+"\t"+x1+"\t"+x2+"\t"+dev+"\t"+dev_min);

							no++;
						}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 1.0 - dev_min / dev_const;

		double[] values = { y0_min, yi_min, f_min, e1_min, s1_min, e2_min, s2_min, r2, dev_const, dev_min };

		if (r2 < R2)
			return null;

		return values;
	}

	public static double[] hillFitFlat(double yl, double yr, double yd, double xl, double xr, double xd,
			boolean[] flags, double[] ws, double[] xs, double[] ys) {
		if (checkRange(flags, ys) == false)
			return null;

		double slope = 1.0;
		double yinf_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 1000000.0;

		for (double yinf = yl; yinf < yr; yinf += yd) {
			for (double x05 = xl; x05 < xr; x05 += xd) {
				double dev = calcHillDeviation(x05, 0.0, yinf, slope, flags, ws, xs, ys);

				if (dev < dev_min) {
					yinf_min = yinf;
					dev_min = dev;
					x05_min = x05;
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);
		double r2 = 1.0 - dev_min / dev_const;

		double[] values = { x05_min, yinf_min, slope, r2, dev_const, dev_min };

		return values;
	}

	public static double calcConstantDeviation(boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double mean = 0.0;

		int no = 0;

		for (int i = 0; i < xs.length; i++)
			if (flags == null || flags[i]) {
				mean += ys[i];
				no++;
			}

		mean /= no;

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++)
			if (flags == null || flags[i]) {
				double delta = (ys[i] - mean);

				if (ws != null)
					dev += ws[i] * delta * delta;
				else
					dev += delta * delta;
			}

		return dev;
	}

	public static double calcHillDeviation(double x05, double y0, double yinf, double slope, boolean[] flags,
			double[] ws, double[] xs, double[] ys) {
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++)

			if (flags == null || flags[i]) {
				double y = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs[i])));
				double delta = (y - ys[i]);

				if (ws != null)
					dev += ws[i] * delta * delta;
				else
					dev += delta * delta;
			}

		return dev;
	}

	public static double calcAUC(double[] fitValues, double x0, double x1, int no) {
		if (fitValues == null)
			return 0.0;

		double auc = 0.0;

		double delta = (x1 - x0) / no;

		double yl = 0.0;

		for (int i = 0; i < no; i++) {
			double xl = x0 + i * delta;
			double xr = x0 + i * delta + delta;

			if (i == 0)
				yl = (fitValues[1] - fitValues[6]) / (1.0 + Math.exp(LN10 * fitValues[2] * (fitValues[0] - xl)));

			double yr = (fitValues[1] - fitValues[6]) / (1.0 + Math.exp(LN10 * fitValues[2] * (fitValues[0] - xr)));

			auc += 0.5 * delta * (yl + yr);

			yl = yr;
		}

		return auc;
	}

	public static double calcAUC(double[] fitValues1, double[] fitValues2, double pct1, double x0, double x1, int no) {
		if (fitValues1 == null || fitValues2 == null)
			return 0.0;

		double auc = 0.0;

		double delta = (x1 - x0) / no;

		double yl = 0.0;

		for (int i = 0; i < no; i++) {
			double xl = x0 + i * delta;
			double xr = x0 + i * delta + delta;

			if (i == 0) {

				double y1 = (fitValues1[1] - fitValues1[6])
						/ (1.0 + Math.exp(LN10 * fitValues1[2] * (fitValues1[0] - xl - Math.log10(pct1))));
				double y2 = (fitValues2[1] - fitValues2[6])
						/ (1.0 + Math.exp(LN10 * fitValues2[2] * (fitValues2[0] - Math.log10(1.0 - pct1) - xl)));
				yl = y1 + y2 - 0.01 * y1 * y2;
			}

			double y1 = (fitValues1[1] - fitValues1[6])
					/ (1.0 + Math.exp(LN10 * fitValues1[2] * (fitValues1[0] - xr - Math.log10(pct1))));
			double y2 = (fitValues2[1] - fitValues2[6])
					/ (1.0 + Math.exp(LN10 * fitValues2[2] * (fitValues2[0] - xr - Math.log10(1 - pct1))));
			double yr = y1 + y2 - 0.01 * y1 * y2;

			auc += 0.5 * delta * (yl + yr);

			yl = yr;
		}

		return auc;
	}

	public static double calcHillY(double[] fitValues, double x) {
		if (x > 0.0)
			x = Math.log10(x);
		else
			return 0.0;

		double y = fitValues[6] + (fitValues[1] - fitValues[6])
				/ (1.0 + Math.exp(LN10 * fitValues[2] * (fitValues[0] - x)));

		return y;
	}

	public static double calcHillY(double x05, double y0, double yinf, double slope, double x) {

		double y = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - x)));

		return y;
	}

	public static double calcHillX(double x05, double y0, double yinf, double slope, double y) {
		double td = (yinf - y0) / (y - y0);
		double x = x05 - Math.log10(td - 1.0) / slope;
		return x;
	}

	public static double calcHillDeviation(double y0, double yinf, double fraction, double ec1, double slope1,
			double ec2, double slope2, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double dev = 0.0;

		double delta1 = (yinf - y0) * fraction;
		double delta2 = (yinf - y0) * (1 - fraction);

		for (int i = 0; i < xs.length; i++)
			if (flags == null || flags[i]) {
				double y = y0 + delta1 / (1.0 + Math.exp(LN10 * slope1 * (ec1 - xs[i]))) + delta2
						/ (1.0 + Math.exp(LN10 * slope2 * (ec2 - xs[i])));

				double delta = (y - ys[i]);

				if (ws != null)
					dev += ws[i] * delta * delta;
				else
					dev += delta * delta;
			}

		return dev;
	}

	public static double[] calcHillFitCurve(double x05, double y0, double yinf, double slope, double[] xs) {
		double[] ys = new double[xs.length];

		for (int i = 0; i < xs.length; i++) {
			ys[i] = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs[i])));
		}

		return ys;
	}

	public static double[] calcHillFitCurve(double y0, double yinf, double fraction, double ec1, double slope1,
			double ec2, double slope2, double[] xs) {
		double delta1 = (yinf - y0) * fraction;
		double delta2 = (yinf - y0) * (1 - fraction);

		double[] ys = new double[xs.length];

		for (int i = 0; i < xs.length; i++) {
			ys[i] = y0 + delta1 / (1.0 + Math.exp(LN10 * slope1 * (ec1 - xs[i]))) + delta2
					/ (1.0 + Math.exp(LN10 * slope2 * (ec2 - xs[i])));
		}

		return ys;
	}

	public static double[] calcHillFitConstant(double[] xs, double[] ys) {
		double[] yy = new double[xs.length];

		double mean = 0.0;
		for (int i = 0; i < xs.length; i++)
			mean += ys[i];
		mean /= xs.length;

		for (int i = 0; i < xs.length; i++) {
			yy[i] = mean;

		}

		return yy;
	}

}