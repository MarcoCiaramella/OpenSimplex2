#define PSIZE 2048
#define PMASK 2047
#define N2 0.01001634121365712
#define N3 0.030485933181293584
#define N4 0.009202377986303158


#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD




typedef struct {
    int xsv, ysv;
    double dx, dy;
} LatticePoint2D;

typedef struct {
    double dx, dy;
} Grad2;

typedef struct {
    short perm[PSIZE];
    Grad2 permGrad2[PSIZE];
} OpenSimplexGradients;

typedef struct {
    Grad2 GRADIENTS_2D[PSIZE];
    LatticePoint2D LOOKUP_2D[4];
} OpenSimplexEnv;



/*
	 * Utility
	 */

int _fastFloor(double x){
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

Grad2 _newGrad2(double dx, double dy){
    Grad2 grad2;
    grad2.dx = dx;
    grad2.dy = dy;
    return grad2;
}

void _loadGrad2ConstArray(OpenSimplexEnv *ose){
    Grad2 arr[24];
    int i = 0;
	arr[i++] = _newGrad2(0.130526192220052, 0.99144486137381);
	arr[i++] = _newGrad2(0.38268343236509, 0.923879532511287);
	arr[i++] = _newGrad2(0.608761429008721, 0.793353340291235);
	arr[i++] = _newGrad2(0.793353340291235, 0.608761429008721);
	arr[i++] = _newGrad2(0.923879532511287, 0.38268343236509);
	arr[i++] = _newGrad2(0.99144486137381, 0.130526192220051);
	arr[i++] = _newGrad2(0.99144486137381, -0.130526192220051);
	arr[i++] = _newGrad2(0.923879532511287, -0.38268343236509);
	arr[i++] = _newGrad2(0.793353340291235, -0.60876142900872);
	arr[i++] = _newGrad2(0.608761429008721, -0.793353340291235);
	arr[i++] = _newGrad2(0.38268343236509, -0.923879532511287);
	arr[i++] = _newGrad2(0.130526192220052, -0.99144486137381);
	arr[i++] = _newGrad2(-0.130526192220052, -0.99144486137381);
	arr[i++] = _newGrad2(-0.38268343236509, -0.923879532511287);
	arr[i++] = _newGrad2(-0.608761429008721, -0.793353340291235);
	arr[i++] = _newGrad2(-0.793353340291235, -0.608761429008721);
	arr[i++] = _newGrad2(-0.923879532511287, -0.38268343236509);
	arr[i++] = _newGrad2(-0.99144486137381, -0.130526192220052);
	arr[i++] = _newGrad2(-0.99144486137381, 0.130526192220051);
	arr[i++] = _newGrad2(-0.923879532511287, 0.38268343236509);
	arr[i++] = _newGrad2(-0.793353340291235, 0.608761429008721);
	arr[i++] = _newGrad2(-0.608761429008721, 0.793353340291235);
	arr[i++] = _newGrad2(-0.38268343236509, 0.923879532511287);
	arr[i++] = _newGrad2(-0.130526192220052, 0.99144486137381);
    for (int i = 0; i < 24; i++) {
        arr[i].dx /= N2;
        arr[i].dy /= N2;
    }
	for (int i = 0; i < PSIZE; i++) {
        ose->GRADIENTS_2D[i] = arr[i % 24];
    }
}

LatticePoint2D _newLatticePoint2D(int xsv, int ysv){
	LatticePoint2D lp2D;
	lp2D.xsv = xsv;
	lp2D.ysv = ysv;
	double ssv = (xsv + ysv) * -0.211324865405187;
	lp2D.dx = -xsv - ssv;
	lp2D.dy = -ysv - ssv;
	return lp2D;
}

void _loadLatticePoint2DConstArray(OpenSimplexEnv *ose){
	ose->LOOKUP_2D[0] = _newLatticePoint2D(1, 0);
	ose->LOOKUP_2D[1] = _newLatticePoint2D(0, 0);
	ose->LOOKUP_2D[2] = _newLatticePoint2D(1, 1);
	ose->LOOKUP_2D[3] = _newLatticePoint2D(0, 1);
}

/*
	 * Noise Evaluators
	 */

/**
	 * 2D Simplex noise base.
	 * Lookup table implementation inspired by DigitalShadow.
	 */
double _noise2_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys){
	double value = 0;

	// Get base points and offsets
	int xsb = _fastFloor(xs), ysb = _fastFloor(ys);
	double xsi = xs - xsb, ysi = ys - ysb;

	// Index to point list
	int index = (int)((ysi - xsi) / 2 + 1);

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi, yi = ysi + ssi;

	// Point contributions
	for (int i = 0; i < 3; i++){
		LatticePoint2D *c = &(ose->LOOKUP_2D[index + i]);

		double dx = xi + c->dx, dy = yi + c->dy;
		double attn = 0.5 - dx * dx - dy * dy;
		if (attn <= 0)
			continue;

		int pxm = (xsb + c->xsv) & PMASK, pym = (ysb + c->ysv) & PMASK;
		Grad2 grad = osg->permGrad2[osg->perm[pxm] ^ pym];
		double extrapolation = grad.dx * dx + grad.dy * dy;

		attn *= attn;
		value += attn * attn * extrapolation;
	}

	return value;
}

/**
	 * 2D Simplex noise, standard lattice orientation.
	 */
double noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y){

	// Get points for A2* lattice
	double s = 0.366025403784439 * (x + y);
	double xs = x + s, ys = y + s;

	return _noise2_Base(ose, osg, xs, ys);
}

/**
	 * 2D Simplex noise, with Y pointing down the main diagonal.
	 * Might be better for a 2D sandbox style game, where Y is vertical.
	 * Probably slightly less optimal for heightmaps or continent maps.
	 */
double noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y){

	// Skew transform and rotation baked into one.
	double xx = x * 0.7071067811865476;
	double yy = y * 1.224744871380249;

	return _noise2_Base(ose, osg, yy + xx, yy - xx);
}

OpenSimplexEnv initOpenSimplex(){
	OpenSimplexEnv ose;
	_loadGrad2ConstArray(&ose);
	_loadLatticePoint2DConstArray(&ose);
	return ose;
}

OpenSimplexGradients newOpenSimplexGradients(OpenSimplexEnv *ose, long seed){
    OpenSimplexGradients osg;
    short source[PSIZE];
    for (short i = 0; i < PSIZE; i++){
        source[i] = i;
    }
    for (int i = PSIZE - 1; i >= 0; i--){
        seed = seed * 6364136223846793005L + 1442695040888963407L;
		int r = (int)((seed + 31) % (i + 1));
		if (r < 0){
            r += (i + 1);
        }
		osg.perm[i] = source[r];
		osg.permGrad2[i] = ose->GRADIENTS_2D[osg.perm[i]];
		source[r] = source[i];
    }
	return osg;
}

__kernel void main(const unsigned int size, __global double* output){
    int x = get_global_id(0);
    int y = get_global_id(1);
    int index = y*get_global_size(0) + x;
    if (index < size){
        OpenSimplexEnv ose = initOpenSimplex();
        OpenSimplexGradients osg = newOpenSimplexGradients(&ose, 1234);
        output[index] = noise2(&ose, &osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ);
    }
}
