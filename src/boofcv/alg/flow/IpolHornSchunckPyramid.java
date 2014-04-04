package boofcv.alg.flow;

// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2012, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.
//
// Java port modifications (C) 2014 by Peter Abeles <peter.abeles@gmail.com>


import boofcv.alg.misc.ImageMiscOps;
import boofcv.struct.image.ImageFloat32;
import java.lang.Math;
import java.lang.System;

/**
 * @author Peter Abeles
 */
public class IpolHornSchunckPyramid
{

	public static final float SOR_EXTRAPOLATION_PARAMETER = 1.9f;
	public static final float INPUT_PRESMOOTHING_SIGMA = 0.8f;

	public static final int DEFAULT_GAUSSIAN_WINDOW_SIZE = 5;

	public static final float ZOOM_SIGMA_ZERO = 0.6f;


	ImageFloat32 I2x = new ImageFloat32(1,1); // x derivative of I2
	ImageFloat32 I2y  = new ImageFloat32(1,1); // y derivative of I2
	ImageFloat32 I2w  = new ImageFloat32(1,1); // warping of I2
	ImageFloat32 I2wx = new ImageFloat32(1,1); // warping of I2x
	ImageFloat32 I2wy = new ImageFloat32(1,1); // warping of I2y
	ImageFloat32 Au   = new ImageFloat32(1,1); // finalant part of numerator of u
	ImageFloat32 Av   = new ImageFloat32(1,1); // finalant part of numerator of v
	ImageFloat32 Du   = new ImageFloat32(1,1); // denominator of u
	ImageFloat32 Dv   = new ImageFloat32(1,1); // denominator of v
	ImageFloat32 D    = new ImageFloat32(1,1); // common numerator of u and v

	/**
	 *
	 *  Function to compute the SOR iteration at a given position
	 *  (SOR = Successive Over-Relaxation)
	 *
	 */
	static float sor_iteration(
			final float []Au, // finalant part of the numerator of u
			final float []Av, // finalant part of the numerator of v
			final float []Du, // denominator of u
			final float []Dv, // denominator of v
			final float []D,  // finalant part of the numerator
			float 	    []u,  // x component of the flow
			float 	    []v,  // y component of the flow
			final float  al, // alpha smoothness parameter
			final int    p,  // current position
			final int    p1, // up-left neighbor
			final int    p2, // up-right neighbor
			final int    p3, // bottom-left neighbor
			final int    p4, // bottom-right neighbor
			final int    p5, // up neighbor
			final int    p6, // left neighbor
			final int    p7, // bottom neighbor
			final int    p8  // right neighbor
	)
	{
		// set the SOR extrapolation parameter
		final float w = SOR_EXTRAPOLATION_PARAMETER;

		// compute the divergence
		final float ula = 1f/12f * (u[p1] + u[p2] + u[p3] + u[p4]) +
			1f/6f  * (u[p5] + u[p6] + u[p7] + u[p8]);
		final float vla = 1f/12f * (v[p1] + v[p2] + v[p3] + v[p4]) +
			1f/6f  * (v[p5] + v[p6] + v[p7] + v[p8]);

		// store the previous values
		final float uk = u[p];
		final float vk = v[p];

		// update the flow
		u[p] = (1.0f - w) * uk + w * (Au[p] - D[p] * v[p] + al * ula) / Du[p];
		v[p] = (1.0f - w) * vk + w * (Av[p] - D[p] * u[p] + al * vla) / Dv[p];

		// return the convergence error
		return (u[p] - uk) * (u[p] - uk) + (v[p] - vk) * (v[p] - vk);
	}



	/**
	 *
	 *  Horn & Schunck method for optical flow estimation at a single scale
	 *
	 */
	public void horn_schunck_optical_flow(
			final ImageFloat32 I1,             // source image
			final ImageFloat32 I2,             // target image
				  ImageFloat32 u,              // x component of optical flow
				  ImageFloat32 v,              // y component of optical flow
			final float  alpha,          // smoothing parameter
			final int    warps,          // number of warpings per scale
			final float  TOL,            // stopping criterion threshold
			final int    maxiter,        // maximum number of iterations
			final boolean verbose         // switch on messages
	)
	{
		int nx = I1.width;
		int ny = I1.height;

		if (verbose) System.err.printf("Single-scale Horn-Schunck of a %dx%d " +
				"image\n\ta=%f nw=%d eps=%f mi=%d v=%s\n", nx, ny,
				alpha, warps, TOL, maxiter, ""+verbose);

		final int   size   = nx * ny;
		final float alpha2 = alpha * alpha;

		//allocate memory
		I2x.reshape(nx,ny);
		I2y.reshape(nx,ny);
		I2w.reshape(nx,ny);
		I2wx.reshape(nx,ny);
		I2wy.reshape(nx,ny);
		Au.reshape(nx,ny);
		Av.reshape(nx,ny);
		Du.reshape(nx,ny);
		Dv.reshape(nx,ny);
		D.reshape(nx,ny);

		// compute the gradient of the second image
		gradient(I2.data, I2x.data, I2y.data,I2.width,I2.height);

		// iterative approximation to the Taylor expansions
		for(int n = 0; n < warps; n++)
		{
			if(verbose) System.err.printf("Warping %d:", n);

			// warp the second image and its derivatives
			IpolBicubicInterpolation.bicubic_interpolation_warp(I2.data, u.data, v.data, I2w.data, nx, ny, true);
			IpolBicubicInterpolation.bicubic_interpolation_warp(I2x.data, u.data, v.data, I2wx.data, nx, ny, true);
			IpolBicubicInterpolation.bicubic_interpolation_warp(I2y.data, u.data, v.data, I2wy.data, nx, ny, true);

			// store the finalant parts of the system
			for(int i = 0; i < size; i++)
			{
				final float I2wl = I2wx.data[i] * u.data[i] + I2wy.data[i] * v.data[i];
				final float dif  = I1.data[i] - I2w.data[i] + I2wl;

				Au.data[i] = dif * I2wx.data[i];
				Av.data[i] = dif * I2wy.data[i];
				Du.data[i] = I2wx.data[i] * I2wx.data[i] + alpha2;
				Dv.data[i] = I2wy.data[i] * I2wy.data[i] + alpha2;
				D.data[i]  = I2wx.data[i] * I2wy.data[i];
			}

			int niter = 0;
			float error = 1000;

			// iterations of the SOR numerical scheme
			while(error > TOL && niter < maxiter)
			{
				niter++;
				error = 0;

				//process the central part of the optical flow
				for(int i = 1; i < ny-1; i++)
					for(int j = 1; j < nx-1; j++)
					{
						final int k = i * nx + j;
						error += sor_iteration(
								Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
								k, k-nx-1, k-nx+1, k+nx-1,
								k+nx+1, k-nx, k-1, k+nx, k+1
						);
					}

				// process the first and last rows
				for(int j = 1; j < nx-1; j++)
				{
					// first row
					int k = j;
					error += sor_iteration(
							Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
							k, k-1, k+1, k+nx-1, k+nx+1,
							k, k-1, k+nx, k+1
					);

					// last row
					k = (ny-1) * nx + j;
					error += sor_iteration(
							Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
							k, k-nx-1, k-nx+1, k-1, k+1,
							k-nx, k-1, k, k+1
					);
				}

				// process the first and last columns
				for(int i = 1; i < ny-1; i++)
				{
					// first column
					int k = i * nx;
					error += sor_iteration(
							Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
							k, k-nx, k-nx+1, k+nx, k+nx+1,
							k-nx, k, k+nx, k+1
					);

					// last column
					k = (i+1) * nx - 1;
					error += sor_iteration(
							Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
							k, k-nx-1, k-nx, k+nx-1, k+nx,
							k-nx, k-1, k+nx, k
					);
				}

				// process the corners
				// up-left corner
				error += sor_iteration(
						Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
						0, 0, 1, nx, nx+1,
						0, 0, nx, 1
				);

				// up-right corner
				int k = nx - 1;
				error += sor_iteration(
						Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
						k, k-1, k, k+nx-1, k+nx,
						k, k-1, k+nx, k
				);

				// bottom-left corner
				k = (ny-1) * nx;
				error += sor_iteration(
						Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
						k, k-nx, k-nx+1,k, k+1,
						k-nx, k, k, k+1
				);

				// bottom-right corner
				k = ny * nx - 1;
				error += sor_iteration(
						Au.data, Av.data, Du.data, Dv.data, D.data, u.data, v.data, alpha2,
						k, k-1, k, k-nx-1, k-nx,
						k-nx, k-1, k, k
				);

				error = (float)Math.sqrt(error / size);
			}

			if(verbose)
				System.err.printf("Iterations %d (%g)\n", niter, error);
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
	 * Function to normalize the images between 0 and 255
	 *
	 **/
	static void image_normalization(
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


	/**
	 *
	 *  Procedure to handle the pyramidal approach.
	 *  This procedure relies on the previous functions to calculate
	 *  large optical flow fields using a pyramidal scheme.
	 *
	 */
	public void horn_schunck_pyramidal(
			final ImageFloat32 I1,              // source image
			final ImageFloat32 I2,              // target image
			ImageFloat32       u,               // x component of optical flow
			ImageFloat32       v,               // y component of optical flow
			final float  alpha,           // smoothing weight
			final int    nscales,         // number of scales
			final float  zfactor,         // zoom factor
			final int    warps,           // number of warpings per scale
			final float  TOL,             // stopping criterion threshold
			final int    maxiter,         // maximum number of iterations
			final boolean  verbose          // switch on messages
	)
	{
		int nx = I1.width;
		int ny = I1.height;

		if (verbose) System.err.printf("Multiscale Horn-Schunck of a %dx%d pair"+
				"\n\ta=%g ns=%d zf=%g nw=%d eps=%g mi=%d\n", nx, ny,
				alpha, nscales, zfactor, warps, TOL, maxiter);

		int size = nx * ny;

		ImageFloat32 []I1s = new ImageFloat32[nscales];
		ImageFloat32 []I2s = new ImageFloat32[nscales];
		ImageFloat32 []us = new ImageFloat32[nscales];
		ImageFloat32 []vs = new ImageFloat32[nscales];


		I1s[0] = new ImageFloat32(nx,ny);
		I2s[0] = new ImageFloat32(nx,ny);

		// normalize the finest scale images between 0 and 255
		image_normalization(I1.data, I2.data, I1s[0].data, I2s[0].data, size);

		// presmoothing the finest scale images
		gaussian(I1s[0].data, nx, ny, INPUT_PRESMOOTHING_SIGMA);
		gaussian(I2s[0].data, nx, ny, INPUT_PRESMOOTHING_SIGMA);

		us[0] = u;
		vs[0] = v;

		// create the scales
		for(int s = 1; s < nscales; s++)
		{
			int nxx_s1 = I1s[s-1].width;
			int nyy_s1 = I1s[s-1].height;

//			zoom_size(nxx_s1, nyy_s1, nxx_s1+s, nyy_s1+s, zfactor);
			int nxx_s = (int)((float) nxx_s1 * zfactor + 0.5);
			int nyy_s = (int)((float) nyy_s1 * zfactor + 0.5);

//			System.out.printf(" pyramid layer size %d %d %d %d %d\n",s,nxx_s1,nyy_s1,nxx_s,nyy_s);

			I1s[s] = new ImageFloat32(nxx_s,nyy_s);
			I2s[s] = new ImageFloat32(nxx_s,nyy_s);
			us[s] = new ImageFloat32(nxx_s,nyy_s);
			vs[s] = new ImageFloat32(nxx_s,nyy_s);

			// compute the zoom from the previous finer scale
			zoom_out(I1s[s-1].data, I1s[s].data, nxx_s1, nyy_s1, zfactor);
			zoom_out(I2s[s-1].data, I2s[s].data, nxx_s1, nyy_s1, zfactor);
		}

		// initialize the flow
		ImageMiscOps.fill(us[nscales-1],0);
		ImageMiscOps.fill(vs[nscales-1],0);

		// pyramidal approximation to the optic flow
		for(int s = nscales-1; s >= 0; s--)
		{
			int nxx_s = I1s[s].width;
			int nyy_s = I1s[s].height;

			if(verbose)
				System.err.printf("Scale: %d %dx%d\n", s, nxx_s, nyy_s);

			// compute the optical flow at this scale
			horn_schunck_optical_flow(
					I1s[s], I2s[s], us[s], vs[s],
					alpha, warps, TOL, maxiter, verbose
			);

			// if this was the last scale, finish now
			if (s==0) break;

			int nxx_s1 = I1s[s-1].width;
			int nyy_s1 = I1s[s-1].height;

			// otherwise, upsample the optical flow

			// zoom the optic flow for the next finer scale
			zoom_in(us[s].data, us[s-1].data, nxx_s, nyy_s, nxx_s1, nyy_s1);
			zoom_in(vs[s].data, vs[s-1].data, nxx_s, nyy_s, nxx_s1, nyy_s1);

			// scale the optic flow with the appropriate zoom factor
			for(int i = 0; i < nxx_s1 * nyy_s1; i++)
			{
				us[s-1].data[i] *= 1.0 / zfactor;
				vs[s-1].data[i] *= 1.0 / zfactor;
			}
		}
	}

	/**
	 *
	 * Compute the gradient of an image using centered differences
	 *
	 */
	void gradient(
			final float []input, // input image
			float []dx,          // computed x derivative
			float []dy,          // computed y derivative
			final int nx,       // image width
			final int ny        // image height
	)
	{
		// compute gradient in the central body of the image
		for(int i = 1; i < ny-1; i++)
		{
			for(int j = 1; j < nx-1; j++)
			{
				final int k = i * nx + j;
				dx[k] = 0.5f*(input[k+1] - input[k-1]);
				dy[k] = 0.5f*(input[k+nx] - input[k-nx]);
			}
		}

		// compute gradient in the first and last rows
		for(int j = 1; j < nx-1; j++)
		{
			dx[j] = 0.5f*(input[j+1] - input[j-1]);
			dy[j] = 0.5f*(input[j+nx] - input[j]);

			final int k = (ny - 1) * nx + j;

			dx[k] = 0.5f*(input[k+1] - input[k-1]);
			dy[k] = 0.5f*(input[k] - input[k-nx]);
		}

		// compute gradient in the first and last columns
		for(int i = 1; i < ny-1; i++)
		{
			final int p = i * nx;
			dx[p] = 0.5f*(input[p+1] - input[p]);
			dy[p] = 0.5f*(input[p+nx] - input[p-nx]);

			final int k = (i+1) * nx - 1;

			dx[k] = 0.5f*(input[k] - input[k-1]);
			dy[k] = 0.5f*(input[k+nx] - input[k-nx]);
		}

		// compute the gradient in the corners
		dx[0] = 0.5f*(input[1] - input[0]);
		dy[0] = 0.5f*(input[nx] - input[0]);

		dx[nx-1] = 0.5f*(input[nx-1] - input[nx-2]);
		dy[nx-1] = 0.5f*(input[2*nx-1] - input[nx-1]);

		dx[(ny-1)*nx] = 0.5f*(input[(ny-1)*nx + 1] - input[(ny-1)*nx]);
		dy[(ny-1)*nx] = 0.5f*(input[(ny-1)*nx] - input[(ny-2)*nx]);

		dx[ny*nx-1] = 0.5f*(input[ny*nx-1] - input[ny*nx-1-1]);
		dy[ny*nx-1] = 0.5f*(input[ny*nx-1] - input[(ny-1)*nx-1]);

	}

	/**
	 *
	 * In-place Gaussian smoothing of an image
	 *
	 */
	void gaussian(
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

	/**
	 *
	 * Downsample an image
	 *
	 **/
	void zoom_out(
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
	void zoom_in(
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

}
