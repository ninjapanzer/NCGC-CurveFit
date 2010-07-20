/*
 * Copyright 2010 The Scripps Research Institute
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

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class NCGCFitFunction implements FitFunction {

	private static double LN10 = Math.log(10D);

	private static FitFunction f = new NCGCFitFunction();
	
	public static FitFunction getInstance() {
		return f;
	}
	
	public double getResponse(Curve curve, double conc) {
		conc = Math.log10(conc);
		double x05 = curve.getLogEC50();
		double y0 = curve.getYZero();
		double yinf = curve.getYInflection();
		double slope = curve.getHillSlope();
		double y = y0 + (yinf - y0) / (1.0D + Math.exp(LN10 * slope * (x05 - conc)));
		return y;
	}
}