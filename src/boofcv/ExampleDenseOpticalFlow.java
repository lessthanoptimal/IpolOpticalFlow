/*
 * Copyright (c) 2011-2014, Peter Abeles. All Rights Reserved.
 *
 * This file is part of BoofCV (http://boofcv.org).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package boofcv;

import boofcv.abst.flow.DenseOpticalFlow;
import boofcv.abst.flow.IpolBroxSpacial_to_DenseOpticalFlow;
import boofcv.abst.flow.IpolHornSchunkPyramid_to_DenseOpticalFlow;
import boofcv.alg.distort.DistortImageOps;
import boofcv.alg.interpolate.TypeInterpolate;
import boofcv.core.image.ConvertBufferedImage;
import boofcv.factory.flow.FactoryDenseOpticalFlow;
import boofcv.gui.PanelGridPanel;
import boofcv.gui.feature.VisualizeOpticalFlow;
import boofcv.gui.image.AnimatePanel;
import boofcv.gui.image.ShowImages;
import boofcv.io.MediaManager;
import boofcv.io.image.UtilImageIO;
import boofcv.io.wrapper.DefaultMediaManager;
import boofcv.struct.flow.ImageFlow;
import boofcv.struct.image.ImageFloat32;
import boofcv.struct.image.ImageUInt8;

import java.awt.image.BufferedImage;
import java.io.FileNotFoundException;

/**
 * TODO Comment
 *
 * @author Peter Abeles
 */
public class ExampleDenseOpticalFlow {

	public static void main(String[] args) throws FileNotFoundException {
		MediaManager media = DefaultMediaManager.INSTANCE;

//		String fileName0 = "images/dogdance07.png";
//		String fileName1 = "images/dogdance08.png";

		String fileName0 = "images/Urban2_07.png";
		String fileName1 = "images/Urban2_08.png";

//		String fileName0 = "images/Grove2_07.png";
//		String fileName1 = "images/Grove2_09.png";

		DenseOpticalFlow<ImageFloat32> denseFlow =
//				new IpolHornSchunkPyramid_to_DenseOpticalFlow();
				new IpolBroxSpacial_to_DenseOpticalFlow();

		BufferedImage buff0 = media.openImage(fileName0);
		BufferedImage buff1 = media.openImage(fileName1);


		// Dense optical flow is very computationally expensive.  Just process the image at 1/2 resolution
		ImageFloat32 previous = new ImageFloat32(buff0.getWidth(),buff0.getHeight());
		ImageFloat32 current = new ImageFloat32(previous.width,previous.height);
		ImageFlow flow = new ImageFlow(previous.width,previous.height);

		ConvertBufferedImage.convertFrom(buff0, previous);
		ConvertBufferedImage.convertFrom(buff1, current);

		// compute dense motion
		long start = System.currentTimeMillis();
		denseFlow.process(previous, current, flow);
		long stop = System.currentTimeMillis();
		System.out.println(" elapsed "+(stop-start));

		UtilOpticalFlow.saveFlow(flow,"denseflow.bflow");

		// Visualize the results
		PanelGridPanel gui = new PanelGridPanel(1,2);

		BufferedImage converted0 = new BufferedImage(current.width,current.height, BufferedImage.TYPE_INT_RGB);
		BufferedImage converted1 = new BufferedImage(current.width,current.height, BufferedImage.TYPE_INT_RGB);
		BufferedImage visualized = new BufferedImage(current.width,current.height, BufferedImage.TYPE_INT_RGB);

		ConvertBufferedImage.convertTo(previous, converted0, true);
		ConvertBufferedImage.convertTo(current, converted1, true);
		VisualizeOpticalFlow.colorized(flow, 10, visualized);

		AnimatePanel animate = new AnimatePanel(150,converted0,converted1);
		gui.add(animate);
		gui.add(visualized);
		animate.start();

		ShowImages.showWindow(gui,"Dense Optical Flow");
	}
}
