/*
 * Copyright 2010 The Scripps Research Institute
 * 
 * This code has been translated from Perl Module: Color::Calc 
 * http://search.cpan.org/~cfaerber/Color-Calc-0.20a/Calc.pm
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
package edu.scripps.fl.commons.color;

import java.awt.Color;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class ColorUtils {

	public static Color invert(Color c) {
		return new Color(255 - c.getRed(), 255 - c.getGreen(), 255 - c.getBlue());
	}

	public static Color toGreyScale(Color c) {
		int rgb = (int) (((double) c.getRed() * 0.299) + ((double) c.getGreen() * 0.587) + ((double) c.getBlue() * 0.114));
		if (rgb > 255)
			rgb = 255;
		return new Color(rgb, rgb, rgb);
	}

	public static Color mix(Color c1, Color c2) {
		return mix(c1, c2, 0.5);
	}

	public static Color mix(Color c1, Color c2, double alpha) {
		return new Color((int) (c1.getRed() + (c2.getRed() - c1.getRed()) * alpha), (int) (c1.getGreen() + (c2.getGreen() - c1.getGreen()) * alpha),
				(int) (c1.getBlue() + (c2.getBlue() - c1.getBlue()) * alpha));

	}

	public static Color contrast(Color c) {
		return contrast(c, 0.5);
	}

	public static Color contrast(Color c, double alpha) {
		if ((alpha < 0.0) || (alpha > 1.0)) {
			alpha = 0.5;
		}

		int contrast = (int) (255 * alpha);

		int[] rgb = new int[] { c.getRed(), c.getGreen(), c.getBlue() };

		for (int i = 0; i < rgb.length; i++) {
			rgb[i] = (rgb[i] > contrast) ? 0 : 255;
		}

		return new Color(rgb[0], rgb[1], rgb[2]);
	}

	public static Color contrastBlackOrWhite(Color c) {
		return contrastBlackOrWhite(c, 0.5);
	}

	public static Color contrastBlackOrWhite(Color c, double alpha) {
		if ((alpha < 0.0) || (alpha > 1.0)) {
			alpha = 0.5;
		}
		return contrast(toGreyScale(c), alpha);
	}

	public static Color blend(Color c) {
		return blend(c, 0.5);
	}

	public static Color blend(Color c, double alpha) {
		return mix(c, contrast(c), alpha);
	}

	public static Color blendBlackOrWhite(Color c) {
		return blendBlackOrWhite(c, 0.5);
	}

	public static Color blendBlackOrWhite(Color c, double alpha) {
		return mix(c, contrastBlackOrWhite(c), alpha);
	}

	public static Color[] getExcelColors() {
		return new Color[] { new Color(255, 0, 0), new Color(0, 255, 0), new Color(0, 0, 255), new Color(255, 255, 0), new Color(255, 0, 255),
				new Color(0, 255, 255), new Color(128, 0, 0), new Color(0, 128, 0), new Color(0, 0, 128), new Color(128, 128, 0),
				new Color(128, 0, 128), new Color(0, 128, 128), new Color(192, 192, 192), new Color(128, 128, 128), new Color(153, 153, 255),
				new Color(153, 51, 102), new Color(255, 255, 204), new Color(204, 255, 255), new Color(102, 0, 102), new Color(255, 128, 128),
				new Color(0, 102, 204), new Color(204, 204, 255), new Color(0, 0, 128), new Color(255, 0, 255), new Color(255, 255, 0),
				new Color(0, 255, 255), new Color(128, 0, 128), new Color(128, 0, 0), new Color(0, 128, 128), new Color(0, 0, 255),
				new Color(0, 204, 255), new Color(204, 255, 255), new Color(204, 255, 204), new Color(255, 255, 153), new Color(153, 204, 255),
				new Color(255, 153, 204), new Color(204, 153, 255), new Color(255, 204, 153), new Color(51, 102, 255), new Color(51, 204, 204),
				new Color(153, 204, 0), new Color(255, 204, 0), new Color(255, 153, 0), new Color(255, 102, 0), new Color(102, 102, 153),
				new Color(150, 150, 150), new Color(0, 51, 102), new Color(51, 153, 102), new Color(0, 51, 0), new Color(51, 51, 0),
				new Color(153, 51, 0), new Color(153, 51, 102), new Color(51, 51, 153), new Color(51, 51, 51) };
	}

	// 0.2 for red <-> yellow part of spectrum.
	public static Color[] getRedtoYellowSpectrum(int length) {
		return getSpectrum(length, 0.0, 0.20);
	}

	public static Color[] getRedtoGreenSpectrum(int length) {
		return getSpectrum(length, 0.0, 0.3);
	}

	// 0.7 for red <-> blue part of spectrum.
	public static Color[] getRedtoBlueSpectrum(int length) {
		return getSpectrum(length, 0.0, 0.7);
	}

	public static Color[] getSpectrum(int length, double lowerBound, double upperBound) {
		if (length < 2)
			length = 2;

		Color[] colors = new Color[length];

		int ii = length;
		double ds = (upperBound - lowerBound) / (double) length;

		for (double dd = 0; dd < (upperBound - lowerBound) && ii > 0; dd += ds) {
			colors[--ii] = Color.getHSBColor((float) dd, 1.0f, 1.0f);
		}
		return colors;
	}

	public static Color[] getColorsBetween(Color c1, Color c2, int length) {
		if (length < 2)
			length = 2;

		Color[] colors = new Color[length];
		colors[0] = c1;
		colors[colors.length - 1] = c2;

		double alpha = 1D / (double) length;

		for (int ii = 1; ii < length - 1; ii++) {
			colors[ii] = mix(c1, c2, alpha * (double) ii);
		}

		return colors;
	}

	public static Color[] getColorsBetween(Color c1, Color c2, Color c3, int length) {
		int len = length / 2;
		Color[] cs1 = getColorsBetween(c1, c2, len);
		Color[] cs2 = getColorsBetween(c2, c3, len);
		Color[] colors = new Color[(cs1.length + cs2.length) - 1];
		System.arraycopy(cs1, 0, colors, 0, cs1.length);
		System.arraycopy(cs2, 1, colors, cs1.length, cs2.length - 1);
		return colors;
	}

	public static String color2Hex(Color c) {
		return Integer.toHexString(c.getRGB() & 0x00ffffff);
	}

	public static Color hex2Color(String hex) {
		return new Color(Integer.valueOf(hex, 16).intValue());
	}

	public static void main(String[] args) {
		Color[] colors = getRedtoBlueSpectrum(10);
		for (int i = 0; i < colors.length; i++) {
			System.out.print(colors[i].getRGB());
		}

		// Color c = Color.GRAY;
		// System.out.println(c);
		// System.out.println( color2Hex(c) );
		// System.out.println( hex2Color( color2Hex(c) ).toString() );

		showTestForm();
	}

	private static void showTestForm() {
		final java.util.List list = new java.util.ArrayList();
		list.add(getColorsBetween(Color.BLACK, Color.WHITE, 40));
		list.add(getColorsBetween(Color.RED, Color.BLUE, 40));
		list.add(getColorsBetween(Color.RED, Color.WHITE, Color.GREEN, 40));
		list.add(getRedtoBlueSpectrum(40));
		list.add(getSpectrum(40, 0.0, 1.0));
		list.add(getRedtoYellowSpectrum(40));
		list.add(getRedtoGreenSpectrum(40));
		list.add(getExcelColors());
		javax.swing.JFrame f = new javax.swing.JFrame();
		f.getContentPane().add(new javax.swing.JPanel() {
			public void paintComponent(java.awt.Graphics g) {
				int w = this.getWidth() / list.size();
				for (int jj = 0; jj < list.size(); jj++) {
					Color[] c = (Color[]) list.get(jj);
					int h = this.getHeight() / c.length;
					for (int ii = 0; ii < c.length; ii++) {
						g.setColor(c[ii]);
						g.fillRect(jj * w, ii * h, w, h);
					}
				}

			}
		});
		f.setDefaultCloseOperation(f.DISPOSE_ON_CLOSE);
		f.setSize(400, 400);
		f.setVisible(true);

	}
}