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

import java.awt.Polygon;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;

import org.apache.commons.beanutils.ConvertUtils;

/**
 * 
 * @author Mark Southern (southern at scripps dot edu)
 * 
 */
public class ShapeFactory {

	class Circle implements ShapeCreator {
		public Shape create(double delta) {
			return new Ellipse2D.Double(-delta, -delta, delta * 2, delta * 2);
		}
	}

	class Diamond implements ShapeCreator {
		public Shape create(double delta) {
			return new Polygon(intArray(0.0, delta, 0.0, -delta), intArray(-delta, 0.0, delta, 0.0), 4);
		}
	}

	class DownTriangle implements ShapeCreator {
		public Shape create(double delta) {
			return new Polygon(intArray(-delta, +delta, 0.0), intArray(-delta, -delta, delta), 3);
		}
	}

	class HorizontalEllipse implements ShapeCreator {
		public Shape create(double delta) {
			return new Ellipse2D.Double(-delta, -delta / 2, delta * 2, delta);
		}
	}

	class HorizontalRectangle implements ShapeCreator {
		public Shape create(double delta) {
			return new Rectangle2D.Double(-delta, -delta / 2, delta * 2, delta);
		}
	}

	class LeftTriangle implements ShapeCreator {
		public Shape create(double delta) {
			return new Polygon(intArray(-delta, delta, delta), intArray(0.0, -delta, +delta), 3);
		}
	}

	class RightTriangle implements ShapeCreator {
		public Shape create(double delta) {
			return new Polygon(intArray(-delta, delta, -delta), intArray(-delta, 0.0, delta), 3);
		}
	}

	interface ShapeCreator {
		public Shape create(double delta);
	}

	class Square implements ShapeCreator {
		public Shape create(double delta) {
			return new Rectangle2D.Double(-delta, -delta, delta * 2, delta * 2);
		}
	}

	class UpTriangle implements ShapeCreator {
		public Shape create(double delta) {
			return new Polygon(intArray(0.0, delta, -delta), intArray(-delta, delta, delta), 3);
		}
	}

	private static int[] intArray(double... values) {
		return (int[]) ConvertUtils.convert(values, int[].class);
	}

	private ShapeCreator[] shapeCreators = new ShapeCreator[] { new Circle(), new UpTriangle(), /*new Diamond(),*/ new HorizontalRectangle(),
			new DownTriangle(), new HorizontalEllipse(), new RightTriangle(), new LeftTriangle(), new Square() };

	private int shapeIndex;

	public Shape getNextShape(double delta) {
		int idx = this.shapeIndex % this.shapeCreators.length;
		Shape result = shapeCreators[idx].create(delta);
		this.shapeIndex++;
		return result;
	}

	public Shape[] getShapes(double delta) {
		Shape[] shapes = new Shape[shapeCreators.length];
		for (int ii = 0; ii < shapes.length; ii++) {
			shapes[ii] = shapeCreators[ii].create(delta);
		}
		return shapes;
	}
}
