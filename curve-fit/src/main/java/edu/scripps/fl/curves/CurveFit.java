/*
 * This class incorporates code from NCGC CurveFit (http://www.ncgc.nih.gov/pub/openhts/curvefit/)
 * a "United States Government Work". Specifically it derives from and makes use of classes within 
 * the gov.nih.ncgc.batch package. 
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.scripps.fl.curves;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.beanutils.ConvertUtils;
import org.apache.commons.collections.CollectionUtils;

import gov.nih.ncgc.batch.BatchHill;
import gov.nih.ncgc.batch.CurveClass;
import gov.nih.ncgc.batch.HillConstants;
import gov.nih.ncgc.batch.HillCurve;
import gov.nih.ncgc.batch.HillFit;
import gov.nih.ncgc.batch.HillStat;

/**
 * @author Mark Southern (southern at scripps dot edu)
 */
public class CurveFit {
	public static double[] expand(double x[]) {
		double x2[] = new double[2 * x.length - 1];
		for (int i = 0; i < x2.length; i++)
			if (i % 2 == 0)
				x2[i] = x[i / 2];
			else
				x2[i] = 0.5D * (x[(i - 1) / 2] + x[(i + 1) / 2]);
		return x2;
	}

	public static void fit(Curve curve) {
		double y[] = (double[]) ConvertUtils.convert(curve.getResponses(), double[].class);
		double x[] = (double[]) ConvertUtils.convert(curve.getConcentrations(), double[].class);
		for (int ii = 0; ii < x.length; ii++)
			x[ii] = Math.log10(x[ii]);
		// max, min and range
		double minY = Double.MAX_VALUE;
		double maxY = -Double.MAX_VALUE;
		double maxResp = y[y.length - 1];
		for (int i = 0; i < y.length; i++) {
			minY = Math.min(minY, y[i]);
			maxY = Math.max(maxY, y[i]);
		}
		curve.setResponseMin(minY);
		curve.setResponseMax(maxY);
		curve.setResponseRange(maxY - minY);
		curve.setMaxResponse(maxResp);
		// fit
		boolean flags[] = null;
		Map maps[] = null;
		double fitValues[] = null;
		Object fitResults[] = HillFit.doHill(x, y, null, HillConstants.FIT_ITER_NO, HillConstants.P4_FIT);
		if (fitResults != null) {
			flags = (boolean[]) fitResults[0];
			fitValues = (double[]) fitResults[1];
			maps = (Map[]) fitResults[2];
			curve.setYZero(fitValues[6]);
			curve.setLogEC50(fitValues[0]);
			curve.setYInflection(fitValues[1]);
			curve.setHillSlope(fitValues[2]);
			curve.setR2(fitValues[3]);

			double ec50 = 1000000D * Math.exp(Math.log(10D) * curve.getLogEC50());
			double testEC50 = Math.pow(10, curve.getLogEC50());
			Double ic50 = null;
			double logIC50 = BatchHill.iccalc(curve.getYZero(), curve.getYInflection(), curve.getLogEC50(), curve.getHillSlope(), 50D);
			if (logIC50 < 0.0D)
				ic50 = 1000000D * Math.exp(Math.log(10D) * logIC50);
			int dn = Math.max(1, x.length - 4);
			double df = dn;
			double p = HillStat.calcPValue(curve.getYZero(), curve.getYInflection(), curve.getLogEC50(), curve.getHillSlope(), x, y, flags);
			int mask = 0;
			for (int i = 0; i < x.length; i++)
				if (!flags[i])
					mask++;
			double ss = HillStat.calcHillDeviation(curve.getLogEC50(), curve.getYZero(), curve.getYInflection(), curve.getHillSlope(),
					flags, null, x, y);
			curve.setEC50(ec50);
			curve.setIC50(ic50);
			curve.setPHill(p);
			curve.setSYX(ss / df);
			for (int ii = 0; ii < flags.length; ii++) {
				if (flags[ii] == true) {
					curve.setMasked(true);
					break;
				}
			}
		} else {
			curve.setLogEC50(null);
			curve.setHillSlope(null);
			curve.setR2(null);
			curve.setYInflection(null);
			curve.setYZero(null);
			curve.setEC50(null);
			curve.setIC50(null);
			curve.setPHill(null);
			curve.setSYX(null);
			curve.setMasked(false);
			flags = new boolean[x.length];
		}
		// masks
		List<Boolean> masks = new ArrayList<Boolean>(flags.length);
		CollectionUtils.addAll(masks, (Boolean[]) ConvertUtils.convert(flags, Boolean[].class));
		curve.setMask(masks);
		// classify
		curveClassification(curve, y, x, flags);
		// rank
		double rank = -BatchHill.calcRank(curve.getCurveClass(), curve.getMaxResponse(), curve.getResponseRange(), fitValues);
		curve.setRank(rank);
	}

	public static void curveClassification(Curve curve) {
		double y[] = (double[]) ConvertUtils.convert(curve.getResponses(), double[].class);
		double x[] = (double[]) ConvertUtils.convert(curve.getConcentrations(), double[].class);
		for (int ii = 0; ii < x.length; ii++)
			x[ii] = Math.log10(x[ii]);
		boolean flags[] = (boolean[]) ConvertUtils.convert(curve.getMask(), boolean[].class);
		curveClassification(curve, y, x, flags);
	}

	protected static void curveClassification(Curve curve, double[] responses, double[] logConcentrations, boolean[] flags) {
		HillCurve hc = new HillCurve(responses, logConcentrations, flags, null);
		if (curve.getHillSlope() != null) // has curve
			hc.setCurve(curve.getLogEC50(), curve.getHillSlope(), curve.getYZero(), curve.getYInflection(), curve.getR2(),
					HillConstants.FIT_ITER_NO);
		CurveClass cc = new CurveClass();
		cc.setSD(HillConstants.CLASSIFICATION_SD);
		cc.setSDfactor(HillConstants.CLASSIFICATION_SD_FACTOR);
		cc.setR2Cutoff(HillConstants.R2);
		cc.setAllowBell(HillConstants.bellMask);
		double curveClass = cc.enzCurveClass(hc);
		curve.setCurveClass(curveClass);
		String signalDirection = getSignalDirection(curveClass);
		curve.setSignalDirection(signalDirection);
		String curveDesc = getCurveDescription(curveClass);
		curve.setCurveDescription(curveDesc);
	}

	protected static String getSignalDirection(double classNo) {
		String signalDirection = "Inactive";
		if (classNo == 4D)
			signalDirection = "Inactive";
		else if (classNo < 0.0D)
			signalDirection = "Decrease";
		else if (classNo > 0.0D)
			signalDirection = "Increase";
		return signalDirection;
	}

	protected static String getCurveDescription(double classNo) {
		String curveDesc = "Inactive";
		if (Math.abs(classNo) == 1.1)
			curveDesc = "Full Curve, >80% POC";
		else if (Math.abs(classNo) == 1.2)
			curveDesc = "Full Curve, <80% POC";
		else if (Math.abs(classNo) == 2.1)
			curveDesc = "Partial Curve, >80% POC";
		else if (Math.abs(classNo) == 2.2)
			curveDesc = "Partial Curve, <80% POC";
		else if (Math.abs(classNo) == 3D)
			curveDesc = "Single Point Active";
		else if (classNo == 4D)
			curveDesc = "Inactive";
		else
			curveDesc = "Noisy Curve";
		return curveDesc;
	}
}