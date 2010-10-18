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
package edu.scripps.fl.curves.plot;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.beanutils.ConvertUtils;
import org.apache.commons.lang.math.NumberUtils;

import com.googlecode.charts4j.AxisLabels;
import com.googlecode.charts4j.AxisLabelsFactory;
import com.googlecode.charts4j.AxisStyle;
import com.googlecode.charts4j.AxisTextAlignment;
import com.googlecode.charts4j.Color;
import com.googlecode.charts4j.Data;
import com.googlecode.charts4j.DataUtil;
import com.googlecode.charts4j.Fills;
import com.googlecode.charts4j.GCharts;
import com.googlecode.charts4j.LegendPosition;
import com.googlecode.charts4j.Plots;
import com.googlecode.charts4j.Shape;
import com.googlecode.charts4j.XYLine;
import com.googlecode.charts4j.XYLineChart;

import edu.scripps.fl.curves.Curve;
import edu.scripps.fl.curves.FitFunction;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class GCurvePlot {
	
	protected static XYLine sampleFunctionToLine(Curve curve, FitFunction f, double start, double end, int samples) {
		double yValues[] = new double[samples];
		double xValues[] = new double[samples];
		
		double step = (end - start) / (double) (samples - 1);
		for (int i = 0; i < samples; i++) {
			double x = start + step * (double) i;
			xValues[i] = x;
			double y = f.getResponse(curve, Math.pow(10, x));
			yValues[i] = y;
		}
		
		Data xData = DataUtil.scaleWithinRange(NumberUtils.min(xValues), NumberUtils.max(xValues), xValues);
		Data yData = DataUtil.scaleWithinRange(NumberUtils.min(yValues), NumberUtils.max(yValues), yValues);
		return Plots.newXYLine(xData, yData, Color.GREEN, "");
	}
	private Color backgroundColor = Color.WHITE;
	private CurvePlotDrawingSupplier drawingSupplier = new CurvePlotDrawingSupplier();
	private List<XYLine> lines = new ArrayList();
	private double minX, maxX, minY, maxY;
	private String title = ""; // "Concentration v Response"
	
	private int width = 500, height = 400;

	public void addCurve(Curve curve, FitFunction function) {
		double[] yValues = (double[]) ConvertUtils.convert(curve.getResponses(), double[].class);
        double curveMinY = NumberUtils.min(yValues);
        double curveMaxY = NumberUtils.max(yValues);
        this.minY = Math.min(minY, curveMinY);
        this.maxY = Math.min(maxY, curveMaxY);
        Data yData = DataUtil.scaleWithinRange(curveMinY, curveMaxY, yValues);
        
        double[] xValues = (double[]) ConvertUtils.convert(curve.getConcentrations(), double[].class);
        for(int ii = 0; ii < xValues.length; ii++) {
        	double x = Math.log10(xValues[ii]);
        	xValues[ii] = x;
        }
        double curveMinX = NumberUtils.min(xValues);
        double curveMaxX = NumberUtils.max(xValues);
        this.minX = Math.min(minX, curveMinX);
        this.maxX = Math.min(maxX, curveMaxX);
        Data xData = DataUtil.scaleWithinRange(NumberUtils.min(xValues), NumberUtils.max(xValues), xValues);
        
        String hexColor = Integer.toHexString( ((java.awt.Color) drawingSupplier.getNextPaint()).getRGB() & 0x00ffffff );
        StringBuffer sb = new StringBuffer();
        sb.append(hexColor);
        while(sb.length() < 6 )
        	sb.insert(0, "0");
        Color color = Color.newColor( sb.toString() );
        
        XYLine line1 = Plots.newXYLine(xData, yData, getBackgroundColor(), "");
//        line1.setLineStyle(LineStyle.newLineStyle(3, 1, 0));
        line1.addShapeMarkers(Shape.CIRCLE, color, 5);
        		
		XYLine fittedLine = sampleFunctionToLine(curve, function, curveMinX, curveMaxX, 100);
//		fittedLine.setLineStyle(LineStyle.newLineStyle(3, 1, 0));
		fittedLine.setColor(color);
		
		lines.add(line1);
		lines.add(fittedLine);
	}

	public Color getBackgroundColor() {
		return backgroundColor;
	}

	public int getHeight() {
		return height;
	}

	public String getURL() {
		XYLineChart chart = GCharts.newXYLineChart(lines);
		chart.setLegendPosition(LegendPosition.BOTTOM);
        chart.setSize(getWidth(), getHeight());
        chart.setTitle(title, Color.BLACK, 14);

        // Defining axis info and styles
        AxisStyle axisStyle = AxisStyle.newAxisStyle(Color.BLACK, 12, AxisTextAlignment.CENTER);
        AxisLabels yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(-25, 125, 25);
        yAxis.setAxisStyle(axisStyle);
        
        AxisLabels xAxis = AxisLabelsFactory.newNumericRangeAxisLabels(Math.floor(minX), Math.ceil(maxX), 1.0);
        xAxis.setAxisStyle(axisStyle);

        chart.addYAxisLabels(yAxis);
        chart.addXAxisLabels(xAxis);
//        chart.setGrid(100, 6.78, 5, 0);

        // Defining background and chart fills.
        chart.setBackgroundFill(Fills.newSolidFill(getBackgroundColor()));
        chart.setAreaFill(Fills.newSolidFill(getBackgroundColor()));
        String url = chart.toURLString();        
        return url;
	}

	public int getWidth() {
		return width;
	}

	public void setBackgroundColor(Color backgroundColor) {
		this.backgroundColor = backgroundColor;
	}
	
	public void setHeight(int height) {
		this.height = height;
	}
	
	public void setWidth(int width) {
		this.width = width;
	}
}