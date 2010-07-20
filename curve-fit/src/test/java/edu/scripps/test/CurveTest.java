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
public class CurveTest {

	public static Curve curve1() {
		return fullCurveIncrease();
	}

	public static Curve fullCurveDecrease() {
		Curve curve = new Curve();
		double[] conc = new double[] { 7.0E-10, 3.7E-9, 1.83E-8, 9.15E-8, 4.572E-7, 2.286E-6, 1.14293E-5 };
		double[] resp = new double[] { 1.9290, 0.5883, -4.8260, -30.0400, -74.8800, -93.0300, -98.3100 };
		for (int ii = 0; ii < conc.length; ii++) {
			curve.add(resp[ii], conc[ii]);
		}
		return curve;
	}

	public static Curve fullCurveIncrease() {
		Curve curve = new Curve();
		double[] conc = new double[] { 3.12E-08, 3.12E-08, 3.12E-08, 6.23E-08, 6.23E-08, 6.23E-08, 1.25E-07, 1.25E-07, 1.25E-07, 2.49E-07,
				2.49E-07, 2.49E-07, 4.98E-07, 4.98E-07, 4.98E-07, 9.97E-07, 9.97E-07, 9.97E-07, 1.99E-06, 1.99E-06, 1.99E-06, 3.99E-06,
				3.99E-06, 3.99E-06, 7.97E-06, 7.97E-06, 7.97E-06, 1.60E-05, 1.60E-05, 1.60E-05 };
		double[] resp = new double[] { 54.274, 6.81, 4.752, 10.905, 13.074, 8.728, 37.15, 28.369, 19.967, 104.683, 102.905, 84.791,
				104.886, 99.521, 70.937, 105.037, 103.962, 105.762, 104.207, 103.253, 105.138, 104.46, 103.125, 105.237, 104.272, 103.359,
				105.287, 103.817, 103.508, 104.777 };
		for (int ii = 0; ii < conc.length; ii++) {
			curve.add(resp[ii], conc[ii]);
		}
		return curve;
	}

	public static Curve inactiveCurve() {
		Curve curve = new Curve();
		double[] conc = new double[] { 3.7E-9, 1.83E-8, 9.15E-8, 4.572E-7, 2.286E-6, 1.14293E-5, 5.71429E-5 };
		double[] resp = new double[] { -0.1437, -0.4878, -0.4988, 0.3642, -0.1635, -0.3200, -0.3996 };
		for (int ii = 0; ii < conc.length; ii++) {
			curve.add(resp[ii], conc[ii]);
		}
		return curve;
	}
	
	public static void main(String[] args) throws Exception {
		testJFreeChartPlots();
		testGooglePlots();
	}

	public static Curve partialCurveDecrease() {
		Curve curve = new Curve();
		double[] conc = new double[] { 3.7E-9, 1.83E-8, 9.15E-8, 4.572E-7, 2.286E-6, 1.14293E-5, 5.71429E-5 };
		double[] resp = new double[] { 0.1501, 0.4061, -0.7567, 1.3190, -4.3990, -29.7100, -66.6800 };
		for (int ii = 0; ii < conc.length; ii++) {
			curve.add(resp[ii], conc[ii]);
		}
		return curve;
	}

	public static void testGooglePlots() throws Exception {
		GCurvePlot gPlot = new GCurvePlot();
		
		Curve curve = CurveTest.fullCurveDecrease();
		CurveFit.fit(curve);
		gPlot.addCurve(curve, new NCGCFitFunction());
		
		curve = CurveTest.partialCurveDecrease();
		CurveFit.fit(curve);
		gPlot.addCurve(curve, new NCGCFitFunction());
		
		String url = gPlot.getURL();
		System.out.println("Generated Google Chart URL: " + url);
	}
	
	public static void testJFreeChartPlots() throws Exception {
		CurvePlot plot = new CurvePlot();
		CurvePlotDrawingSupplier d = (CurvePlotDrawingSupplier) plot.getDrawingSupplier();
		d.setLineWidth(1);
		d.setShapeSize(2);
		plot.setResponseRange(-120, 20);
		for (Curve curve : new Curve[] { fullCurveDecrease(), partialCurveDecrease() } ) {
//		for (Curve curve : new Curve[] { inactiveCurve(), fullCurveDecrease(), partialCurveDecrease(), fullCurveIncrease() }) {
			CurveFit.fit(curve);
			System.out.println(curve);
			plot.addCurveMeanAndStdDev(curve, new NCGCFitFunction());
		}
		plot.addLineAt(-50);
		File file = File.createTempFile("curveplot", ".png");
		System.out.println("Writing file " + file);
		plot.write(new FileOutputStream(file));
		Desktop.getDesktop().open(file);
	}
	
	public static void testSerialization(Curve curve) throws IOException {
		File file = File.createTempFile("aaa", "xml");
		file.deleteOnExit();
		java.beans.XMLEncoder enc = new XMLEncoder(new FileOutputStream(file));
		enc.writeObject(curve);
		enc.close();
		java.beans.XMLDecoder dec = new XMLDecoder(new FileInputStream(file));
		curve = (Curve) dec.readObject();
		dec.close();
	}
}