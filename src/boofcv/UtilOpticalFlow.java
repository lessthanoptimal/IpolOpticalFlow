package boofcv;

import boofcv.struct.flow.ImageFlow;

import java.io.*;

/**
 * @author Peter Abeles
 */
public class UtilOpticalFlow {

	public static void saveFlow( ImageFlow image , String fileName )
			throws FileNotFoundException
	{
		PrintStream out = new PrintStream(fileName);

		out.println(image.width);
		out.println(image.height);
		for (int y = 0; y < image.height; y++) {
			for (int x = 0; x < image.width; x++) {
				ImageFlow.D d = image.unsafe_get(x,y);

				out.print(d.x+" "+d.y+" ");
			}
			out.println();
		}
		out.close();
	}

	public static ImageFlow loadFlow( String fileName )
			throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(fileName));

		int width = Integer.parseInt(in.readLine());
		int height = Integer.parseInt(in.readLine());

		ImageFlow image = new ImageFlow(width,height);

		for (int y = 0; y < image.height; y++) {
			String line = in.readLine();
			String words[] = line.split("\\s");
			for (int x = 0; x < image.width; x++) {
				ImageFlow.D d = image.unsafe_get(x,y);

				d.x = Float.parseFloat(words[x*2]);
				d.y = Float.parseFloat(words[x*2+1]);
			}
		}

		return image;
	}
}
