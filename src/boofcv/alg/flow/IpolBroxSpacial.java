package boofcv.alg.flow;

import boofcv.alg.InputSanityCheck;
import boofcv.struct.image.ImageFloat32;

import java.util.List;

import static boofcv.alg.flow.IpolBicubicInterpolation.bicubic_interpolation_warp;
import static boofcv.alg.flow.IpolDerivative.divergence_u;
import static boofcv.alg.flow.IpolDerivative.psi_divergence;

/**
 * @author Peter Abeles
 */
public class IpolBroxSpacial
{
	public static final double EPSILON = 0.001;
	public static final double MAXITER = 300;
	public static final double SOR_PARAMETER = 1.9;
	public static final double GAUSSIAN_SIGMA = 0.8;

	//allocate memory
	float[] du    = new float[1];
	float[] dv    = new float[1];

	float[] ux    = new float[1];
	float[] uy    = new float[1];
	float[] vx    = new float[1];
	float[] vy    = new float[1];

	float[] I1x   = new float[1];
	float[] I1y   = new float[1];
	float[] I2x   = new float[1];
	float[] I2y   = new float[1];
	float[] I2w   = new float[1];
	float[] I2wx  = new float[1];
	float[] I2wy  = new float[1];
	float[] I2xx  = new float[1];
	float[] I2yy  = new float[1];
	float[] I2xy  = new float[1];
	float[] I2wxx = new float[1];
	float[] I2wyy = new float[1];
	float[] I2wxy = new float[1];

	float[] div_u = new float[1];
	float[] div_v = new float[1];
	float[] div_d = new float[1];

	float[] Au    = new float[1];
	float[] Av    = new float[1];
	float[] Du    = new float[1];
	float[] Dv    = new float[1];
	float[] D     = new float[1];

	float[] psid  = new float[1];
	float[] psig  = new float[1];
	float[] psis  = new float[1];
	float[] psi1  = new float[1];
	float[] psi2  = new float[1];
	float[] psi3  = new float[1];
	float[] psi4  = new float[1];

	/**
	 *
	 * Compute the coefficients of the robust functional (data term)
	 *
	 **/
	public static void psi_data(
			final float []I1,  //first image
			final float []I2,  //second image
			final float []I2x, //gradient of the second image
			final float []I2y, //gradient of the second image
			final float []du,  //motion increment
			final float []dv,  //motion increment
			float []psip,      //output coefficients
			final int nx,     //image width
			final int ny      //image height
	)
	{
		final int size = nx * ny;

		//compute 1/(sqrt((I2-I1+I2x*du+I2y*dv)²+e²) in each pixel
		for(int i = 0; i < size; i++)
		{
			final float dI  = I2[i] - I1[i] + I2x[i] * du[i] + I2y[i] * dv[i];
			final float dI2 = dI * dI;

			psip[i] = (float)(1.0 / Math.sqrt(dI2 + EPSILON * EPSILON));
		}
	}

	/**
	 *
	 * Compute the coefficients of the robust functional (gradient term)
	 *
	 **/
	public static void psi_gradient(
			final float []I1x,  //gradient of the first image
			final float []I1y,  //gradient of the first image
			final float []I2x,  //gradient of the second image
			final float []I2y,  //gradient of the second image
			final float []I2xx, //second derivatives of the second image
			final float []I2xy, //second derivatives of the second image
			final float []I2yy, //second derivatives of the second image
			final float []du,   //motion increment
			final float []dv,   //motion increment
			float []psip,       //output coefficients
			final int nx,      //image width
			final int ny       //image height
	)
	{
		final int size = nx * ny;

		//compute 1/(sqrt(|DI2-DI1+HI2*(du,dv)|²+e²) in each pixel
		for(int i = 0; i < size; i++)
		{
			final float dIx = I2x[i] - I1x[i] + I2xx[i] * du[i] + I2xy[i] * dv[i];
			final float dIy = I2y[i] - I1y[i] + I2xy[i] * du[i] + I2yy[i] * dv[i];
			final float dI2 = dIx * dIx + dIy * dIy;

			psip[i] = (float)(1. / Math.sqrt(dI2 + EPSILON * EPSILON));
		}
	}

	/**
	 *
	 * Compute the coefficients of the robust functional (smoothness term)
	 *
	 **/
	public static void psi_smooth(
			final float []ux, //gradient of x component of the optical flow
			final float []uy, //gradient of x component of the optical flow
			final float []vx, //gradient of y component of the optical flow
			final float []vy, //gradient of y component of the optical flow
			float []psi,      //output coefficients
			final int nx,    //image width
			final int ny     //image height
	)
	{
		final int size = nx * ny;

		//compute 1/(sqrt(ux²+uy²+vx²+vy²+e²) in each pixel
		for(int i = 0; i < size; i++)
		{
			final float du  = ux[i] * ux[i] + uy[i] * uy[i];
			final float dv  = vx[i] * vx[i] + vy[i] * vy[i];
			final float d2  = du + dv;

			psi[i] = (float)(1. / Math.sqrt(d2 + EPSILON * EPSILON));
		}
	}

	/**
	 *
	 *  SOR iteration in one position
	 *
	 */
	public static float sor_iteration(
			final float []Au,   //finalant part of the numerator of u
			final float []Av,   //finalant part of the numerator of v
			final float []Du,   //denominator of u
			final float []Dv,   //denominator of v
			final float []D,    //finalant part of the numerator
			float       []du,   //x component of the motion increment
			float       []dv,   //y component of the motion increment
			final float alpha, //alpha smoothness parameter
			final float []psi1, //coefficients of the divergence
			final float []psi2,
			final float []psi3,
			final float []psi4,
			final int   i,     //current row
			final int   i0,    //previous row
			final int   i1,    //following row
			final int   j,     //current column
			final int   nx,    //number of columns
			final int   j0,    //previous column
			final int   j1     //following column
	)
	{
		//set the SOR extrapolation parameter
		final float w = (float)SOR_PARAMETER;

		//calculate the position in the array
		final int k = i * nx + j;

		//compute the divergence part of the numerator
		final float div_du = psi1[k] * du[k+i1] + psi2[k] * du[k-i0] +
				psi3[k] * du[k+j1] + psi4[k] * du[k-j0] ;
		final float div_dv = psi1[k] * dv[k+i1] + psi2[k] * dv[k-i0] +
				psi3[k] * dv[k+j1] + psi4[k] * dv[k-j0] ;

		final float duk = du[k];
		final float dvk = dv[k];

		//update the motion increment
		du[k] = (1.0f-w) * du[k] + w * (Au[k] - D[k] * dv[k] + alpha * div_du) / Du[k];
		dv[k] = (1.0f-w) * dv[k] + w * (Av[k] - D[k] * du[k] + alpha * div_dv) / Dv[k];

		//return the covergence error in this position
		return (du[k] - duk) * (du[k] - duk) + (dv[k] - dvk) * (dv[k] - dvk);
	}

	public void declareMemory( int w, int h ) {
		int N = w*h;
		if( du.length >= N )
			return;

		du = new float[N];
		dv = new float[N];

		ux = new float[N];
		uy = new float[N];
		vx = new float[N];
		vy = new float[N];

		I1x = new float[N];
		I1y = new float[N];
		I2x = new float[N];
		I2y = new float[N];
		I2w = new float[N];
		I2wx = new float[N];
		I2wy = new float[N];
		I2xx = new float[N];
		I2yy = new float[N];
		I2xy = new float[N];
		I2wxx = new float[N];
		I2wyy = new float[N];
		I2wxy = new float[N];

		div_u = new float[N];
		div_v = new float[N];
		div_d = new float[N];

		Au = new float[N];
		Av = new float[N];
		Du = new float[N];
		Dv = new float[N];
		D = new float[N];

		psid = new float[N];
		psig = new float[N];
		psis = new float[N];
		psi1 = new float[N];
		psi2 = new float[N];
		psi3 = new float[N];
		psi4 = new float[N];
	}

	/**
	 *
	 * Compute the optic flow with the Brox spatial method
	 *
	 **/
	void brox_optic_flow
	(
			final float[] I1,   //first image
			final float[] I2,   //second image
			float[] u,          //x component of the optical flow
			float[] v,          //y component of the optical flow
			final int    nx,         //image width
			final int    ny,         //image height
			final float  alpha,      //smoothness parameter
			final float  gamma,      //gradient term parameter
			final float  TOL,        //stopping criterion threshold
			final int    inner_iter, //number of inner iterations
			final int    outer_iter, //number of outer iterations
			final boolean   verbose  //switch on messages
	)
	{
		final int size = nx * ny;

		declareMemory(nx,ny);

		//compute the gradient of the images
		IpolDerivative.gradient(I1, I1x, I1y, nx, ny);
		IpolDerivative.gradient(I2, I2x, I2y, nx, ny);

		//compute second order derivatives
		IpolDerivative.Dxx(I2, I2xx, nx, ny);
		IpolDerivative.Dyy(I2, I2yy, nx, ny);
		IpolDerivative.Dxy(I2, I2xy, nx, ny);

		//outer iterations loop
		for(int no = 0; no < outer_iter; no++)
		{
			//warp the second image and its derivatives
			bicubic_interpolation_warp(I2, u, v, I2w, nx, ny, true);
			bicubic_interpolation_warp(I2x, u, v, I2wx, nx, ny, true);
			bicubic_interpolation_warp(I2y, u, v, I2wy, nx, ny, true);
			bicubic_interpolation_warp(I2xx, u, v, I2wxx, nx, ny, true);
			bicubic_interpolation_warp(I2xy, u, v, I2wxy, nx, ny, true);
			bicubic_interpolation_warp(I2yy, u, v, I2wyy, nx, ny, true);

			//compute the flow gradient
			IpolDerivative.gradient(u, ux, uy, nx, ny);
			IpolDerivative.gradient(v, vx, vy, nx, ny);

			//compute robust function Phi for the smoothness term
			psi_smooth(ux, uy, vx, vy, psis, nx, ny);

			//compute coefficients of Phi functions in divergence
			psi_divergence(psis, psi1, psi2, psi3, psi4, nx, ny);

			//compute the divergence for the gradient of w
			divergence_u(u, v, psi1, psi2, psi3, psi4, div_u, div_v, nx, ny);

			for(int i = 0; i < size; i++)
			{
				//compute the coefficents of dw[i] in the smoothness term
				div_d[i] = alpha * (psi1[i] + psi2[i] + psi3[i] + psi4[i]);

				//initialize the motion increment
				du[i] = dv[i] = 0;
			}

			//inner iterations loop
			for(int ni = 0; ni < inner_iter; ni++)
			{
				//compute robust function Phi for the data and gradient terms
				psi_data(I1, I2w, I2wx, I2wy, du, dv,  psid, nx, ny);
				psi_gradient(I1x, I1y, I2wx, I2wy, I2wxx, I2wxy, I2wyy, du, dv, psig, nx, ny);

				//store finalant parts of the numerical scheme
				for(int i = 0; i < size; i++)
				{
					final float p = psid[i];
					final float g = gamma * psig[i];

					//brightness finalancy term
					final float dif = I2w[i] - I1[i];
					final float BNu = -p * dif * I2wx[i];
					final float BNv = -p * dif * I2wy[i];
					final float BDu = p * I2wx[i] * I2wx[i];
					final float BDv = p * I2wy[i] * I2wy[i];

					//gradient finalancy term
					final float dx  = (I2wx[i] - I1x[i]);
					final float dy  = (I2wy[i] - I1y[i]);
					final float GNu = -g * (dx * I2wxx[i] + dy * I2wxy[i]);
					final float GNv = -g * (dx * I2wxy[i] + dy * I2wyy[i]);
					final float GDu =  g * (I2wxx[i] * I2wxx[i] + I2wxy[i] * I2wxy[i]);
					final float GDv =  g * (I2wyy[i] * I2wyy[i] + I2wxy[i] * I2wxy[i]);
					final float DI  = (I2wxx[i] + I2wyy[i]) * I2wxy[i];
					final float Duv =  p * I2wy[i] * I2wx[i] + g * DI;

					Au[i] = BNu + GNu + alpha * div_u[i];
					Av[i] = BNv + GNv + alpha * div_v[i];
					Du[i] = BDu + GDu + div_d[i];
					Dv[i] = BDv + GDv + div_d[i];
					D [i] = Duv;
				}

				//sor iterations loop
				float error = 1000;
				int nsor = 0;

				while( error > TOL && nsor < MAXITER)
				{
					error = 0;
					nsor++;

					//update the motion increment in the center of the images
					for(int i = 1; i < ny-1; i++)
						for(int j = 1; j < nx-1; j++) {

							error += sor_iteration(
									Au, Av, Du, Dv, D, du, dv, alpha,
									psi1, psi2, psi3, psi4,
									i, nx, nx, j, nx, 1, 1
							);
						}

					//update the motion increment in the first and last rows
					for(int j = 1; j < nx-1; j++)
					{
						error += sor_iteration(
								Au, Av, Du, Dv, D, du, dv, alpha,
								psi1, psi2, psi3, psi4,
								0, 0, nx, j, nx, 1, 1
						);

						error += sor_iteration(
								Au, Av, Du, Dv, D, du, dv, alpha,
								psi1, psi2, psi3, psi4,
								ny-1, nx, 0, j, nx, 1, 1
						);
					}

					//update the motion increment in the first and last columns
					for(int i = 1; i < ny-1; i++)
					{
						error += sor_iteration(
								Au, Av, Du, Dv, D, du, dv, alpha,
								psi1, psi2, psi3, psi4,
								i, nx, nx, 0, nx, 0, 1
						);

						error += sor_iteration(
								Au, Av, Du, Dv, D, du, dv, alpha,
								psi1, psi2, psi3, psi4,
								i, nx, nx, nx-1, nx, 1, 0
						);
					}

					//process the top-left corner (0,0)
					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							0, 0, nx, 0, nx, 0, 1
					);

					//process the top-right corner (0,nx-1)
					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							0, 0, nx, nx-1, nx, 1, 0
					);

					//process the bottom-left corner (ny-1,0)
					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							ny-1, nx, 0, 0, nx, 0, 1
					);

					//process the bottom-right corner (ny-1,nx-1)
					error += sor_iteration(
							Au, Av, Du, Dv, D, du, dv, alpha,
							psi1, psi2, psi3, psi4,
							ny-1, nx, 0, nx-1, nx, 1, 0
					);

					error = (float)Math.sqrt(error / size);
				}

				if(verbose) System.out.println("Iterations: " + nsor );
			}

			//update the flow with the estimated motion increment
			for(int i = 0; i < size; i++)
			{
				u[i] += du[i];
				v[i] += dv[i];
			}
		}
	}

	/**
	 *
	 *  Multiscale approach for computing the optical flow
	 *
	 **/
	public void brox_optic_flow(
			final ImageFloat32 I1,   //first image
			final ImageFloat32 I2,   //second image
			ImageFloat32 u, 		    //x component of the optical flow
			ImageFloat32 v, 		    //y component of the optical flow
			final float  alpha,      //smoothness parameter
			final float  gamma,      //gradient term parameter
			final int    nscales,    //number of scales
			final float  nu,         //downsampling factor
			final float  TOL,        //stopping criterion threshold
			final int    inner_iter, //number of inner iterations
			final int    outer_iter, //number of outer iterations
			final boolean   verbose  //switch on messages
	)
	{
		int nxx = I1.width;
		int nyy = I1.height;
		int size = nxx * nyy;

		float[][] I1s = new float[nscales][];
		float[][] I2s = new float[nscales][];
		float[][] us = new float[nscales][];
		float[][] vs = new float[nscales][];

		int[] nx = new int[nscales];
		int[] ny = new int[nscales];

		I1s[0] = new float[size];
		I2s[0] = new float[size];

		//normalize the input images between 0 and 255
		IpolFlowUtil.image_normalization(I1.data, I2.data, I1s[0], I2s[0], size);

		//presmoothing the finest scale images
		IpolFlowUtil.gaussian(I1s[0], nxx, nyy, GAUSSIAN_SIGMA);
		IpolFlowUtil.gaussian(I2s[0], nxx, nyy, GAUSSIAN_SIGMA);

		us [0] = u.data;
		vs [0] = v.data;
		nx [0] = nxx;
		ny [0] = nyy;

		//create the scales
		for(int s = 1; s < nscales; s++)
		{
			nx[s] = (int)((float) nx[s-1] * nu + 0.5f);
			ny[s] = (int)((float) ny[s-1] * nu + 0.5f);
			final int sizes = nx[s] * ny[s];

			I1s[s] = new float[sizes];
			I2s[s] = new float[sizes];
			us[s]  = new float[sizes];
			vs[s]  = new float[sizes];

			//compute the zoom from the previous scale
			IpolFlowUtil.zoom_out(I1s[s-1], I1s[s], nx[s-1], ny[s-1], nu);
			IpolFlowUtil.zoom_out(I2s[s-1], I2s[s], nx[s-1], ny[s-1], nu);
		}

		//initialization of the optical flow at the coarsest scale
		for(int i = 0; i < nx[nscales-1] * ny[nscales-1]; i++)
			us[nscales-1][i] = vs[nscales-1][i] = 0.0f;


		//pyramidal approach for computing the optical flow
		for(int s = nscales-1; s >= 0; s--)
		{
			if(verbose) System.out.println("Scale: "+s);

			//compute the optical flow for the current scale
			brox_optic_flow(
					I1s[s], I2s[s], us[s], vs[s], nx[s], ny[s],
					alpha, gamma, TOL, inner_iter, outer_iter, verbose
			);

			//if it is not the finer scale, then upsample the optical flow and adapt it conveniently
			if(s != 0)
			{
				IpolFlowUtil.zoom_in(us[s], us[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);
				IpolFlowUtil.zoom_in(vs[s], vs[s-1], nx[s], ny[s], nx[s-1], ny[s-1]);

				for(int i = 0; i < nx[s-1] * ny[s-1]; i++)
				{
					us[s-1][i] *= 1.0 / nu;
					vs[s-1][i] *= 1.0 / nu;
				}
			}
		}
	}

}
