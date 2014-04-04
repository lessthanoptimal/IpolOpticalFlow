package boofcv.abst.flow;

import boofcv.alg.InputSanityCheck;
import boofcv.alg.flow.IpolHornSchunckPyramid;
import boofcv.struct.flow.ImageFlow;
import boofcv.struct.image.ImageFloat32;
import boofcv.struct.image.ImageType;

/**
 *
 * #define PAR_DEFAULT_NPROC 0
 * #define PAR_DEFAULT_ALPHA 7
 * #define PAR_DEFAULT_NSCALES 10
 * #define PAR_DEFAULT_ZFACTOR 0.5
 * #define PAR_DEFAULT_NWARPS 10
 * #define PAR_DEFAULT_TOL 0.0001
 * #define PAR_DEFAULT_MAXITER 150
 * #define PAR_DEFAULT_VERBOSE 0
 * #define PAR_MAX_ZFACTOR 0.99
 *
 * @author Peter Abeles
 */
public class IpolHornSchunkPyramid_to_DenseOpticalFlow implements DenseOpticalFlow<ImageFloat32> {

	IpolHornSchunckPyramid alg = new IpolHornSchunckPyramid();

	float alpha;
	int nscales;
	float zfactor;
	int warps;
	float TOL;
	int maxIter;

	ImageFloat32 u = new ImageFloat32(1,1);
	ImageFloat32 v = new ImageFloat32(1,1);

	/**
	 *
	 * @param alpha Try 7.
	 * @param nscales Number of scales in image pyramid.  try 10
	 * @param zfactor Try 0.5
	 * @param warps Try 10
	 * @param TOL Tolerance for convergence
	 * @param maxIter
	 */
	public IpolHornSchunkPyramid_to_DenseOpticalFlow(float alpha, int nscales, float zfactor, int warps, float TOL, int maxIter) {
		this.alpha = alpha;
		this.nscales = nscales;
		this.zfactor = zfactor;
		this.warps = warps;
		this.TOL = TOL;
		this.maxIter = maxIter;
	}

	public IpolHornSchunkPyramid_to_DenseOpticalFlow() {
		this.alpha = 7;
		this.nscales = 10;
		this.zfactor = 0.5f;
		this.warps = 10;
		this.TOL = 0.0001f;
		this.maxIter = 150;
	}

	@Override
	public void process(ImageFloat32 source, ImageFloat32 destination, ImageFlow flow) {

		InputSanityCheck.checkSameShape(source,destination);

		int nx = source.width;
		int ny = source.height;

		// Set the number of scales according to the size of the
		// images.  The value N is computed to assure that the smaller
		// images of the pyramid don't have a size smaller than 16x16
		int nscales = this.nscales;
		final int N = (int)(1 + Math.log(Math.hypot(nx, ny)/16) / Math.log(1/zfactor));
		if(N < nscales)
			nscales = N;

		// it can't handle sub-images
		if( source.isSubimage() )
			source = source.clone();
		if( destination.isSubimage() )
			destination = destination.clone();

		u.reshape(source.width,source.height);
		v.reshape(source.width,source.height);

		alg.horn_schunck_pyramidal(source,destination,u,v,alpha,nscales,zfactor,warps,TOL,maxIter,false);

		for( int y = 0; y < flow.height; y++ ) {
			for( int x = 0; x < flow.width; x++ ) {
				ImageFlow.D f = flow.unsafe_get(x,y);
				f.x = u.unsafe_get(x,y);
				f.y = v.unsafe_get(x,y);
			}
		}
	}

	@Override
	public ImageType<ImageFloat32> getInputType() {
		return ImageType.single(ImageFloat32.class);
	}
}
