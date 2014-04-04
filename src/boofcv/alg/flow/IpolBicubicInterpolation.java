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

/**
 * @author Peter Abeles
 */
public class IpolBicubicInterpolation {
	public static int BOUNDARY_CONDITION = 0;
//0 Neumann
//1 PerioBoundaryTypedic
//2 Symmetric

	/**
	 *
	 * Neumann boundary condition test
	 *
	 **/
	static int neumann_bc(int x, int nx, BoundaryHit out)
	{
		if(x < 0)
		{
			x = 0;
			out.hit = true;
		}
		else if (x >= nx)
		{
			x = nx - 1;
			out.hit = true;
		}

		return x;
	}

	/**
	 *
	 * Periodic boundary condition test
	 *
	 **/
	static int periodic_bc(int x, int nx, BoundaryHit out)
	{
		if(x < 0)
		{
			final int n   = 1 - (int)(x/(nx+1));
			final int ixx = x + n * nx;

			x =   ixx% nx;
			out.hit = true;
		}
		else if(x >= nx)
		{
			x = x % nx;
			out.hit = true;
		}

		return x;
	}


	/**
	 *
	 * Symmetric boundary condition test
	 *
	 **/
	static int symmetric_bc(int x, int nx, BoundaryHit out)
	{
		if(x < 0)
		{
			final int borde = nx - 1;
			final int xx = -x;
			final int n  = (int)(xx/borde) % 2;

			if ( n != 0 ) x = borde - ( xx % borde );
			else x = xx % borde;
			out.hit = true;
		}
		else if ( x >= nx )
		{
			final int borde = nx - 1;
			final int n = (int)(x/borde) % 2;

			if ( n != 0) x = borde - ( x % borde );
			else x = x % borde;
			out.hit = true;
		}

		return x;
	}


	/**
	 *
	 * Cubic interpolation in one dimension
	 *
	 **/
	static double cubic_interpolation_cell (
			double v[],  //interpolation points.  size 4
			double x      //point to be interpolated
	)
	{
		return  v[1] + 0.5 * x * (v[2] - v[0] +
				x * (2.0 *  v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] +
						x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
	}


	/**
	 *
	 * Bicubic interpolation in two dimensions
	 *
	 **/
	static double bicubic_interpolation_cell (
			double p[][], //array 4x4 containing the interpolation points
			double x,       //x position to be interpolated
			double y        //y position to be interpolated
	)
	{
		double v[] = new double[4];
		v[0] = cubic_interpolation_cell(p[0], y);
		v[1] = cubic_interpolation_cell(p[1], y);
		v[2] = cubic_interpolation_cell(p[2], y);
		v[3] = cubic_interpolation_cell(p[3], y);
		return cubic_interpolation_cell(v, x);
	}

	/**
	 *
	 * Compute the bicubic interpolation of a point in an image.
	 * Detect if the point goes outside the image domain.
	 *
	 **/
	public static float bicubic_interpolation_at(
			final float []input, //image to be interpolated
			final float  uu,    //x component of the vector field
			final float  vv,    //y component of the vector field
			final int    nx,    //image width
			final int    ny,    //image height
			boolean       border_out //if true, return zero outside the region
	)
	{
		final int sx = (uu < 0)? -1: 1;
		final int sy = (vv < 0)? -1: 1;

		int x, y, mx, my, dx, dy, ddx, ddy;
		BoundaryHit out = new BoundaryHit();
		out.hit = false;

		//apply the corresponding boundary conditions
		switch(BOUNDARY_CONDITION) {

			case 0: x   = neumann_bc((int) uu, nx, out);
				y   = neumann_bc((int) vv, ny, out);
				mx  = neumann_bc((int) uu - sx, nx, out);
				my  = neumann_bc((int) vv - sx, ny, out);
				dx  = neumann_bc((int) uu + sx, nx, out);
				dy  = neumann_bc((int) vv + sy, ny, out);
				ddx = neumann_bc((int) uu + 2*sx, nx, out);
				ddy = neumann_bc((int) vv + 2*sy, ny, out);
				break;

			case 1: x   = periodic_bc((int) uu, nx, out);
				y   = periodic_bc((int) vv, ny, out);
				mx  = periodic_bc((int) uu - sx, nx, out);
				my  = periodic_bc((int) vv - sx, ny, out);
				dx  = periodic_bc((int) uu + sx, nx, out);
				dy  = periodic_bc((int) vv + sy, ny, out);
				ddx = periodic_bc((int) uu + 2*sx, nx, out);
				ddy = periodic_bc((int) vv + 2*sy, ny, out);
				break;

			case 2: x   = symmetric_bc((int) uu, nx, out);
				y   = symmetric_bc((int) vv, ny, out);
				mx  = symmetric_bc((int) uu - sx, nx, out);
				my  = symmetric_bc((int) vv - sx, ny, out);
				dx  = symmetric_bc((int) uu + sx, nx, out);
				dy  = symmetric_bc((int) vv + sy, ny, out);
				ddx = symmetric_bc((int) uu + 2*sx, nx, out);
				ddy = symmetric_bc((int) vv + 2*sy, ny, out);
				break;

			default:x   = neumann_bc((int) uu, nx, out);
				y   = neumann_bc((int) vv, ny, out);
				mx  = neumann_bc((int) uu - sx, nx, out);
				my  = neumann_bc((int) vv - sx, ny, out);
				dx  = neumann_bc((int) uu + sx, nx, out);
				dy  = neumann_bc((int) vv + sy, ny, out);
				ddx = neumann_bc((int) uu + 2*sx, nx, out);
				ddy = neumann_bc((int) vv + 2*sy, ny, out);
				break;
		}

		if( out.hit && border_out)
			return 0.0f;

		else
		{
			//obtain the interpolation points of the image
			final float p11 = input[mx  + nx * my];
			final float p12 = input[x   + nx * my];
			final float p13 = input[dx  + nx * my];
			final float p14 = input[ddx + nx * my];

			final float p21 = input[mx  + nx * y];
			final float p22 = input[x   + nx * y];
			final float p23 = input[dx  + nx * y];
			final float p24 = input[ddx + nx * y];

			final float p31 = input[mx  + nx * dy];
			final float p32 = input[x   + nx * dy];
			final float p33 = input[dx  + nx * dy];
			final float p34 = input[ddx + nx * dy];

			final float p41 = input[mx  + nx * ddy];
			final float p42 = input[x   + nx * ddy];
			final float p43 = input[dx  + nx * ddy];
			final float p44 = input[ddx + nx * ddy];

			//create array
			double pol[][] = new double[][]{
			{p11, p21, p31, p41},
			{p12, p22, p32, p42},
			{p13, p23, p33, p43},
			{p14, p24, p34, p44}
		};

			//return interpolation
			return (float)bicubic_interpolation_cell(pol, uu-x, vv-y);
		}
	}


	/**
	 *
	 * Compute the bicubic interpolation of an image.
	 *
	 **/
	public static void bicubic_interpolation_warp(
			final float []input,  //image to be warped
			final float []u,      //x component of the vector field
			final float []v,      //y component of the vector field
			float       []output, //warped output image with bicubic interpolation
			final int    nx,     //image width
			final int    ny,     //image height
			boolean      border_out//if true, put zeros outside the region
	)
	{
		for(int i = 0; i < ny; i++)
			for(int j = 0; j < nx; j++)
			{
				final int   p  = i * nx + j;
				final float uu = (float) (j + u[p]);
				final float vv = (float) (i + v[p]);

				//obtain the bicubic interpolation at position (uu, vv)
				output[p] = bicubic_interpolation_at(input,
						uu, vv, nx, ny, border_out);
			}
	}
	
	public static class BoundaryHit {
		boolean hit;
	}

}
