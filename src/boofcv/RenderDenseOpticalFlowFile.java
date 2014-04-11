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

import boofcv.gui.feature.VisualizeOpticalFlow;
import boofcv.gui.image.ShowImages;
import boofcv.struct.flow.ImageFlow;

import java.awt.image.BufferedImage;
import java.io.IOException;

/**
 * TODO Comment
 *
 * @author Peter Abeles
 */
public class RenderDenseOpticalFlowFile {

	public static void main(String[] args) throws IOException {

		ImageFlow flow = UtilOpticalFlow.loadFlow("denseflow.bflow");

		// Visualize the results
		BufferedImage visualized = new BufferedImage(flow.width,flow.height, BufferedImage.TYPE_INT_RGB);

		VisualizeOpticalFlow.colorized(flow, 10, visualized);

		ShowImages.showWindow(visualized,"Dense Optical Flow");
	}
}
