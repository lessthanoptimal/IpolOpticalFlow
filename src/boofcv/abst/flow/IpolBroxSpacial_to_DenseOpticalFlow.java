package boofcv.abst.flow;

import boofcv.alg.InputSanityCheck;
import boofcv.alg.flow.IpolBroxSpacial;
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
public class IpolBroxSpacial_to_DenseOpticalFlow implements DenseOpticalFlow<ImageFloat32> {

	IpolBroxSpacial alg = new IpolBroxSpacial();

	float  alpha;      //smoothness parameter
	float  gamma;      //gradient term parameter
	int    nscales;    //number of scales
	float  nu;         //downsampling factor
	float  TOL;        //stopping criterion threshold
	int    inner_iter; //number of inner iterations
	int    outer_iter; //number of outer iterations
	float zfactor;     // reduction factor for creating the scales

	ImageFloat32 u = new ImageFloat32(1,1);
	ImageFloat32 v = new ImageFloat32(1,1);

	public IpolBroxSpacial_to_DenseOpticalFlow(float alpha, int gamma, int nscales , float nu , float zfactor , float TOL , int inner_iter , int outer_iter ) {
		this.alpha = alpha;
		this.gamma = gamma;
		this.nscales = nscales;
		this.nu = nu;
		this.zfactor = zfactor;
		this.TOL = TOL;
		this.inner_iter = inner_iter;
		this.outer_iter = outer_iter;
	}

	public IpolBroxSpacial_to_DenseOpticalFlow() {
		this.alpha = 18;
		this.gamma = 7;
		this.nscales = 100;
		this.nu = 0.75f;
		this.zfactor = 0.75f;
		this.TOL = 0.0001f;
		this.inner_iter = 1;
		this.outer_iter = 15;
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

		alg.brox_optic_flow(source, destination, u, v, alpha, gamma, nscales, nu, TOL, inner_iter, outer_iter, false);

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
