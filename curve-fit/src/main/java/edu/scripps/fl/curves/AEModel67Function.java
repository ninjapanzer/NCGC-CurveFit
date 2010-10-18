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
public class AEModel67Function implements FitFunction {

	private static FitFunction f = new AEModel67Function();
	
	public static FitFunction getInstance() {
		return f;
	}
	public double getResponse(Curve curve, double conc) {
		double Bmax = curve.getYInflection();
		double n = curve.getHillSlope();
		double K = curve.getEC50();
		double Y2 = curve.getYZero();
		return ((Bmax * Math.pow(conc, n)) / (Math.pow(K, n) + Math.pow(conc, n))) + Y2;
	}
}