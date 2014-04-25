package boofcv.alg.flow;

/**
 * @author Peter Abeles
 */
public class IpolFlowUtil
{
	public static final int DEFAULT_GAUSSIAN_WINDOW_SIZE = 5;

	public static final float ZOOM_SIGMA_ZERO = 0.6f;

	/**
	 *
	 * Function to normalize the images between 0 and 255
	 *
	 **/
	public static void image_normalization(
			final float []I1,   // first image
			final float []I2,   // second image
			float       []I1n,  // first normalized image
			float       []I2n,  // second normalized image
			int          size  // size of each image
	)
	{
		// find the max and min of both images
		final float max1 = max_element(I1, size);
		final float max2 = max_element(I2, size);
		final float min1 = min_element(I1, size);
		final float min2 = min_element(I2, size);

		// obtain the absolute max and min
		final float max = max1 > max2 ? max1 : max2;
		final float min = min1 < min2 ? min1 : min2;
		final float den = max - min;

		if(den > 0)
			// normalize both images
			for(int i = 0; i < size; i++)
			{
				I1n[i] = 255.0f * (I1[i] - min) / den;
				I2n[i] = 255.0f * (I2[i] - min) / den;
			}

		else
			// copy the original images
			for(int i = 0; i < size; i++)
			{
				I1n[i] = I1[i];
				I2n[i] = I2[i];
			}
	}

	// compute the largest number of an array
	static float max_element(final float []x, int n)
	{
		int r = 0;
		for (int i = 1; i < n; i++)
			if (x[i] > x[r])
				r = i;
		return x[r];
	}

	// compute the smallest number of an array
	static float min_element(final float []x, int n)
	{
		int r = 0;
		for (int i = 1; i < n; i++)
			if (x[i] < x[r])
				r = i;
		return x[r];
	}

	/**
	 *
	 * Downsample an image
	 *
	 **/
	public static void zoom_out(
			final float []I,          // input image
			float []Iout,             // output image
			final int nx,            // image width
			final int ny,            // image height
			final float factor       // zoom factor between 0 and 1
	)
	{
		// temporary working image
		float []Is = new float[nx*ny];
		System.arraycopy(I,0,Is,0,nx*ny);

		// compute the size of the zoomed image
		int nxx = (int)((float) nx * factor + 0.5);
		int nyy = (int)((float) ny * factor + 0.5);

		// compute the Gaussian sigma for smoothing
		final float sigma = ZOOM_SIGMA_ZERO * (float)Math.sqrt(1.0/(factor*factor) - 1.0);

		// pre-smooth the image
		gaussian(Is, nx, ny, sigma);

		// re-sample the image using bicubic interpolation
		for (int i1 = 0; i1 < nyy; i1++)
			for (int j1 = 0; j1 < nxx; j1++)
			{
				final float i2  = (float) i1 / factor;
				final float j2  = (float) j1 / factor;

				float g = IpolBicubicInterpolation.bicubic_interpolation_at(Is, j2, i2, nx, ny, false);
				Iout[i1 * nxx + j1] = g;
			}
	}


	/**
	 *
	 * Function to upsample the image
	 *
	 **/
	public static void zoom_in(
			final float []I, // input image
			float []Iout,    // output image
			int nx,         // width of the original image
			int ny,         // height of the original image
			int nxx,        // width of the zoomed image
			int nyy         // height of the zoomed image
	)
	{
		// compute the zoom factor
		final float factorx = ((float)nxx / nx);
		final float factory = ((float)nyy / ny);

		// re-sample the image using bicubic interpolation
		for (int i1 = 0; i1 < nyy; i1++)
			for (int j1 = 0; j1 < nxx; j1++)
			{
				float i2 =  (float) i1 / factory;
				float j2 =  (float) j1 / factorx;

				float g = IpolBicubicInterpolation.bicubic_interpolation_at(I, j2, i2, nx, ny, false);
				try{
					Iout[i1 * nxx + j1] = g;
				} catch ( RuntimeException e ) {
					System.out.println();
				}
			}
	}

	/**
	 *
	 * In-place Gaussian smoothing of an image
	 *
	 */
	public static void gaussian(
			float I[],             // input/output image
			final int xdim,       // image width
			final int ydim,       // image height
			final double sigma    // Gaussian sigma
	)
	{
		final BoundaryType boundary_condition = BoundaryType.DEFAULT;
		final int window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

		final double den  = 2*sigma*sigma;
		final int   size = (int) (window_size * sigma) + 1 ;
		final int   bdx  = xdim + size;
		final int   bdy  = ydim + size;

		if (boundary_condition != null && size > xdim) {
			System.err.printf("GaussianSmooth: sigma too large\n");
			System.exit(0);
		}

		// compute the coefficients of the 1D convolution kernel
		double B[] = new double[size];
		for(int i = 0; i < size; i++)
			B[i] = 1 / (sigma * Math.sqrt(2.0 * 3.1415926)) * Math.exp(-i * i / den);

		// normalize the 1D convolution kernel
		double norm = 0;
		for(int i = 0; i < size; i++)
			norm += B[i];
		norm *= 2;
		norm -= B[0];
		for(int i = 0; i < size; i++)
			B[i] /= norm;

		// convolution of each line of the input image
		double R[] = new double[size + xdim + size];

		for (int k = 0; k < ydim; k++)
		{
			int i, j;
			for (i = size; i < bdx; i++)
				R[i] = I[k * xdim + i - size];

			switch (boundary_condition)
			{
				case DIRICHLET:
					for(i = 0, j = bdx; i < size; i++, j++)
						R[i] = R[j] = 0;
					break;

				case REFLECTING:
					for(i = 0, j = bdx; i < size; i++, j++) {
						R[i] = I[k * xdim + size-i];
						R[j] = I[k * xdim + xdim-i-1];
					}
					break;

				case PERIODIC:
					for(i = 0, j = bdx; i < size; i++, j++) {
						R[i] = I[k * xdim + xdim-size+i];
						R[j] = I[k * xdim + i];
					}
					break;
			}

			for (i = size; i < bdx; i++)
			{
				double sum = B[0] * R[i];
				for (j = 1; j < size; j++ )
					sum += B[j] * ( R[i-j] + R[i+j] );
				I[k * xdim + i - size] = (float)sum;
			}
		}

		// convolution of each column of the input image
		double T[] = new double[ size + ydim + size ];

		for (int k = 0; k < xdim; k++)
		{
			int i, j;
			for (i = size; i < bdy; i++)
				T[i] = I[(i - size) * xdim + k];

			switch (boundary_condition)
			{
				case DIRICHLET:
					for (i = 0, j = bdy; i < size; i++, j++)
						T[i] = T[j] = 0;
					break;

				case REFLECTING:
					for (i = 0, j = bdy; i < size; i++, j++) {
						T[i] = I[(size-i) * xdim + k];
						T[j] = I[(ydim-i-1) * xdim + k];
					}
					break;

				case PERIODIC:
					for( i = 0, j = bdx; i < size; i++, j++) {
						T[i] = I[(ydim-size+i) * xdim + k];
						T[j] = I[i * xdim + k];
					}
					break;
			}

			for (i = size; i < bdy; i++)
			{
				double sum = B[0] * T[i];
				for (j = 1; j < size; j++ )
					sum += B[j] * (T[i-j] + T[i+j]);
				I[(i - size) * xdim + k] = (float)sum;
			}
		}
	}
}
