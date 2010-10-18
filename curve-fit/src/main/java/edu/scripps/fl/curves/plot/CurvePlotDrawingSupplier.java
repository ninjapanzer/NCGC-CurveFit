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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jfree.chart.ChartColor;
import org.jfree.chart.plot.DrawingSupplier;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class CurvePlotDrawingSupplier implements DrawingSupplier {

	private int outlinePaintIndex;
	private transient Paint[] outlinePaintSequence = new Paint[] { Color.lightGray };

	private int outlineStrokeIndex;
	private transient Stroke[] outlineStrokeSequence = new Stroke[] { new BasicStroke(1.0f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_BEVEL) };;

	private int paintIndex;
	private Paint[] paintSequence = null; // initialized below

	private ShapeFactory shapeFactory = new ShapeFactory();
	private int shapeIndex;

	private double shapeSize = 12;

	private int strokeIndex;
	private transient Stroke[] strokeSequence = new Stroke[] { new BasicStroke(1.0f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_BEVEL) };
	
	private boolean grayScale = false;

	public boolean isGrayScale() {
		return grayScale;
	}

	public void setGrayScale(boolean grayScale) {
		this.grayScale = grayScale;
	}

	public CurvePlotDrawingSupplier() {
		List<Paint> paints = new ArrayList();
		paints.add(Color.GREEN);
		paints.add(Color.BLACK);
		paints.addAll(Arrays.asList(ChartColor.createDefaultPaintArray()));
		paintSequence = paints.toArray(new Paint[0]);

	}

	@Override
	public Paint getNextFillPaint() {
		Paint result = this.paintSequence[this.paintIndex % this.paintSequence.length];
		this.paintIndex++;
		return result;

	}

	@Override
	public Paint getNextOutlinePaint() {
		Paint result = this.outlinePaintSequence[this.outlinePaintIndex % this.outlinePaintSequence.length];
		this.outlinePaintIndex++;
		return result;

	}

	@Override
	public Stroke getNextOutlineStroke() {
		Stroke result = this.outlineStrokeSequence[this.outlineStrokeIndex % this.outlineStrokeSequence.length];
		this.outlineStrokeIndex++;
		return result;

	}

	@Override
	public Paint getNextPaint() {
		Paint result = this.paintSequence[this.paintIndex % this.paintSequence.length];
		this.paintIndex++;
		if( isGrayScale() )
			result = toGreyScale((Color) result);
		return result;

	}
	
	private static Color toGreyScale(Color c) {
		int rgb = (int) (((double) c.getRed() * 0.299) + ((double) c.getGreen() * 0.587) + ((double) c.getBlue() * 0.114));
		if (rgb > 255)
			rgb = 255;
		return new Color(rgb, rgb, rgb);
	}

	@Override
	public Shape getNextShape() {
		return shapeFactory.getNextShape(getShapeSize() / 2);
	}

	@Override
	public Stroke getNextStroke() {
		Stroke result = this.strokeSequence[this.strokeIndex % this.strokeSequence.length];
		this.strokeIndex++;
		return result;

	}

	public double getShapeSize() {
		return shapeSize;
	}

	public void setShapeSize(double shapeSize) {
		this.shapeSize = shapeSize;
	}
	
	public void setLineWidth(double width) {
		this.strokeSequence = new Stroke[] {
				new BasicStroke((float)width, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_BEVEL)
		};
	}
}