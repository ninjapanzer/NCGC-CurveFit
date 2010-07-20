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

import org.jfree.chart.renderer.xy.*;

import java.awt.*;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.*;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.*;
import org.jfree.data.Range;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.io.SerialUtilities;
import org.jfree.util.BooleanList;
import org.jfree.util.BooleanUtilities;
import org.jfree.util.ObjectUtilities;
import org.jfree.util.PaintUtilities;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class MyXYErrorRenderer extends XYLineAndShapeRenderer
{

    public MyXYErrorRenderer()
    {
        super(false, true);
        errorPaint = null;
        errorStroke = null;
        capLength = 4D;
    }
    
    public boolean isBaseSeriesYError()
    {
		return baseSeriesYError;
	}

	public void setBaseSeriesYError(boolean baseSeriesYError)
	{
		this.baseSeriesYError = baseSeriesYError;
	}

	public boolean isBaseSeriesXError()
	{
		return baseSeriesXError;
	}

	public void setBaseSeriesXError(boolean baseSeriesXError)
	{
		this.baseSeriesXError = baseSeriesXError;
	}
    
    public boolean getSeriesYError(int series)
    {
    	Boolean yError = seriesYError.getBoolean(series);
    	if( yError != null )
    		return yError.booleanValue();
    	return baseSeriesYError;
    }
    
    public void setSeriesYError(int series, boolean flag)
    {
        setSeriesYError(series, BooleanUtilities.valueOf(flag));
    }

    public void setSeriesYError(int series, Boolean flag)
    {
        seriesYError.setBoolean(series, flag);
        fireChangeEvent();
    }
    
    public boolean getSeriesXError(int series) {
    	Boolean xError = seriesXError.getBoolean(series);
    	if( xError != null )
    		return xError.booleanValue();
    	return baseSeriesXError;
    }
    
    public void setSeriesXError(int series, boolean flag)
    {
        setSeriesXError(series, BooleanUtilities.valueOf(flag));
    }

    public void setSeriesXError(int series, Boolean flag)
    {
        seriesXError.setBoolean(series, flag);
        fireChangeEvent();
    }

    public double getCapLength()
    {
        return capLength;
    }

    public void setCapLength(double length)
    {
        capLength = length;
        fireChangeEvent();
    }

    public Paint getErrorPaint()
    {
        return errorPaint;
    }

    public void setErrorPaint(Paint paint)
    {
        errorPaint = paint;
        fireChangeEvent();
    }

    public Stroke getErrorStroke()
    {
        return errorStroke;
    }

    public void setErrorStroke(Stroke stroke)
    {
        errorStroke = stroke;
        fireChangeEvent();
    }

    public Range findDomainBounds(XYDataset dataset)
    {
        if(dataset != null)
            return DatasetUtilities.findDomainBounds(dataset, true);
        else
            return null;
    }

    public Range findRangeBounds(XYDataset dataset)
    {
        if(dataset != null)
            return DatasetUtilities.findRangeBounds(dataset, true);
        else
            return null;
    }

    public void drawItem(Graphics2D g2, XYItemRendererState state, Rectangle2D dataArea, PlotRenderingInfo info, XYPlot plot, ValueAxis domainAxis, ValueAxis rangeAxis, 
            XYDataset dataset, int series, int item, CrosshairState crosshairState, int pass)
    {
        if(pass == 0 && (dataset instanceof IntervalXYDataset) && getItemVisible(series, item))
        {
            IntervalXYDataset ixyd = (IntervalXYDataset)dataset;
            PlotOrientation orientation = plot.getOrientation();
            if(getSeriesXError(series))
            {
                double x0 = ixyd.getStartXValue(series, item);
                double x1 = ixyd.getEndXValue(series, item);
                double y = ixyd.getYValue(series, item);
                org.jfree.ui.RectangleEdge edge = plot.getDomainAxisEdge();
                double xx0 = domainAxis.valueToJava2D(x0, dataArea, edge);
                double xx1 = domainAxis.valueToJava2D(x1, dataArea, edge);
                double yy = rangeAxis.valueToJava2D(y, dataArea, plot.getRangeAxisEdge());
                Line2D cap1 = null;
                Line2D cap2 = null;
                double adj = capLength / 2D;
                Line2D line;
                if(orientation == PlotOrientation.VERTICAL)
                {
                    line = new java.awt.geom.Line2D.Double(xx0, yy, xx1, yy);
                    cap1 = new java.awt.geom.Line2D.Double(xx0, yy - adj, xx0, yy + adj);
                    cap2 = new java.awt.geom.Line2D.Double(xx1, yy - adj, xx1, yy + adj);
                } else
                {
                    line = new java.awt.geom.Line2D.Double(yy, xx0, yy, xx1);
                    cap1 = new java.awt.geom.Line2D.Double(yy - adj, xx0, yy + adj, xx0);
                    cap2 = new java.awt.geom.Line2D.Double(yy - adj, xx1, yy + adj, xx1);
                }
                if(errorPaint != null)
                    g2.setPaint(errorPaint);
                else
                    g2.setPaint(getItemPaint(series, item));
                if(errorStroke != null)
                    g2.setStroke(errorStroke);
                else
                    g2.setStroke(getItemStroke(series, item));
                g2.draw(line);
                g2.draw(cap1);
                g2.draw(cap2);
            }
            if(getSeriesYError(series))
            {
                double y0 = ixyd.getStartYValue(series, item);
                double y1 = ixyd.getEndYValue(series, item);
                double x = ixyd.getXValue(series, item);
                org.jfree.ui.RectangleEdge edge = plot.getRangeAxisEdge();
                double yy0 = rangeAxis.valueToJava2D(y0, dataArea, edge);
                double yy1 = rangeAxis.valueToJava2D(y1, dataArea, edge);
                double xx = domainAxis.valueToJava2D(x, dataArea, plot.getDomainAxisEdge());
                Line2D cap1 = null;
                Line2D cap2 = null;
                double adj = capLength / 2D;
                Line2D line;
                if(orientation == PlotOrientation.VERTICAL)
                {
                    line = new java.awt.geom.Line2D.Double(xx, yy0, xx, yy1);
                    cap1 = new java.awt.geom.Line2D.Double(xx - adj, yy0, xx + adj, yy0);
                    cap2 = new java.awt.geom.Line2D.Double(xx - adj, yy1, xx + adj, yy1);
                } else
                {
                    line = new java.awt.geom.Line2D.Double(yy0, xx, yy1, xx);
                    cap1 = new java.awt.geom.Line2D.Double(yy0, xx - adj, yy0, xx + adj);
                    cap2 = new java.awt.geom.Line2D.Double(yy1, xx - adj, yy1, xx + adj);
                }
                if(errorPaint != null)
                    g2.setPaint(errorPaint);
                else
                    g2.setPaint(getItemPaint(series, item));
                if(errorStroke != null)
                    g2.setStroke(errorStroke);
                else
                    g2.setStroke(getItemStroke(series, item));
                g2.draw(line);
                g2.draw(cap1);
                g2.draw(cap2);
            }
        }
        super.drawItem(g2, state, dataArea, info, plot, domainAxis, rangeAxis, dataset, series, item, crosshairState, pass);
    }

    public boolean equals(Object obj)
    {
        if(obj == this)
            return true;
        if(!(obj instanceof MyXYErrorRenderer))
            return false;
        MyXYErrorRenderer that = (MyXYErrorRenderer)obj;
        if(capLength != that.capLength)
            return false;
        if(!PaintUtilities.equal(errorPaint, that.errorPaint))
            return false;
        if(!ObjectUtilities.equal(errorStroke, that.errorStroke))
            return false;
        else
            return super.equals(obj);
    }

    private void readObject(ObjectInputStream stream)
        throws IOException, ClassNotFoundException
    {
        stream.defaultReadObject();
        errorPaint = SerialUtilities.readPaint(stream);
        errorStroke = SerialUtilities.readStroke(stream);
    }

    private void writeObject(ObjectOutputStream stream)
        throws IOException
    {
        stream.defaultWriteObject();
        SerialUtilities.writePaint(errorPaint, stream);
        SerialUtilities.writeStroke(errorStroke, stream);
    }

    static final long serialVersionUID = 5162283570955172424L;
    private double capLength;
    private transient Paint errorPaint;
    private transient Stroke errorStroke;
    private BooleanList seriesYError = new BooleanList();
    private BooleanList seriesXError = new BooleanList();
    private boolean baseSeriesYError = true;
    private boolean baseSeriesXError = true;
}