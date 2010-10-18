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
package edu.scripps.test;

import java.awt.Desktop;
import java.beans.XMLDecoder;
import java.beans.XMLEncoder;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;

import edu.scripps.fl.curves.AEModel67Function;
import edu.scripps.fl.curves.Curve;
import edu.scripps.fl.curves.CurveFit;
import edu.scripps.fl.curves.NCGCFitFunction;
import edu.scripps.fl.curves.plot.CurvePlot;
import edu.scripps.fl.curves.plot.CurvePlotDrawingSupplier;
import edu.scripps.fl.curves.plot.GCurvePlot;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class TestCurve {

	

	public static Curve curve() {
		Curve curve = new Curve();
		double[] conc = new double[] { 0.00000000297, 0.00000008200, 0.00000036500, 0.00000190600, 0.00001004000, 0.00005361000 };
		double[] resp = new double[] { -117.255, -111.495, -99.4224, -97.8107, -99.9145, -102.267 };
		for (int ii = 0; ii < conc.length; ii++) {
			curve.add(resp[ii], conc[ii]);
		}
		return curve;
	}

	public static void main(String[] args) throws Exception {
		Curve curve = TestCurve.curve();
		CurveFit.fit(curve);
		System.out.println(curve);
	}
}