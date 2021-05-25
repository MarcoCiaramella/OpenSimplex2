#define PMASK 2047
#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD





typedef struct {
    unsigned int length;
    void *data;
} vect;

typedef struct {
    int xsv, ysv;
	double dx, dy;
} LatticePoint2D;

typedef struct {
    double dxr, dyr, dzr;
	int xrv, yrv, zrv;
    int nextOnFailure;
    int nextOnSuccess;
} LatticePoint3D;

typedef struct {
    int xsv, ysv, zsv, wsv;
	double dx, dy, dz, dw;
} LatticePoint4D;

typedef struct {
    double dx, dy;
} Grad2;

typedef struct {
    double dx, dy, dz;
} Grad3;

typedef struct {
    double dx, dy, dz, dw;
} Grad4;




int get_index(const unsigned int width, const unsigned int height){
	int x = get_global_id(0);
	int y = get_global_id(1);
	if (x < width && y < height){
		return y*width + x;
	}
	return -1;
}

double get_x(){
	int x = get_global_id(0);
	return (x + OFF_X) * FREQ;
}

double get_y(){
	int y = get_global_id(1);
	return (y + OFF_Y) * FREQ;
}

double get_z(){
	return 0.0;
}

double get_w(){
	return 0.0;
}

int fast_floor(double x){
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

/**
	 * 2D SuperSimplex noise base.
	 * Lookup table implementation inspired by DigitalShadow.
	 */
double _noise2_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys){
	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs), ysb = fast_floor(ys);
	double xsi = xs - xsb, ysi = ys - ysb;

	// Index to point list
	int a = (int)(xsi + ysi);
	int index =
		(a << 2) |
		(int)(xsi - ysi / 2 + 1 - a / 2.0) << 3 |
		(int)(ysi - xsi / 2 + 1 - a / 2.0) << 4;

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi, yi = ysi + ssi;

	// Point contributions
	for (int i = 0; i < 4; i++){
		LatticePoint2D *c = ose->LOOKUP_2D[index + i];

		double dx = xi + c->dx, dy = yi + c->dy;
		double attn = 2.0 / 3.0 - dx * dx - dy * dy;
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
	 * 2D SuperSimplex noise, standard lattice orientation.
	 */
double noise2(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y){

	// Get points for A2* lattice
	double s = 0.366025403784439 * (x + y);
	double xs = x + s, ys = y + s;

	return _noise2_Base(ose, osg, xs, ys);
}

/**
	 * 2D SuperSimplex noise, with Y pointing down the main diagonal.
	 * Might be better for a 2D sandbox style game, where Y is vertical.
	 * Probably slightly less optimal for heightmaps or continent maps.
	 */
double noise2_XBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y){

	// Skew transform and rotation baked into one.
	double xx = x * 0.7071067811865476;
	double yy = y * 1.224744871380249;

	return _noise2_Base(ose, osg, yy + xx, yy - xx);
}

/**
	 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
	 * Lookup table implementation inspired by DigitalShadow.
	 * It was actually faster to narrow down the points in the loop itself,
	 * than to build up the index with enough info to isolate 8 points.
	 */
double _noise3_BCC(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xr, double yr, double zr){

	// Get base and offsets inside cube of first lattice.
	int xrb = fast_floor(xr), yrb = fast_floor(yr), zrb = fast_floor(zr);
	double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	LatticePoint3D *c = ose->LOOKUP_3D[index];
	while (c != NULL){
		double dxr = xri + c->dxr, dyr = yri + c->dyr, dzr = zri + c->dzr;
		double attn = 0.75 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0){
			c = c->nextOnFailure;
		}
		else{
			int pxm = (xrb + c->xrv) & PMASK, pym = (yrb + c->yrv) & PMASK, pzm = (zrb + c->zrv) & PMASK;
			Grad3 grad = osg->permGrad3[osg->perm[osg->perm[pxm] ^ pym] ^ pzm];
			double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

			attn *= attn;
			value += attn * attn * extrapolation;
			c = c->nextOnSuccess;
		}
	}
	return value;
}

/**
	 * 3D Re-oriented 8-point BCC noise, classic orientation
	 * Proper substitute for what 3D SuperSimplex would be,
	 * in light of Forbidden Formulae.
	 * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
	 */
double noise3_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z){

	// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
	// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
	// Orthonormal rotation. Not a skew transform.
	double r = (2.0 / 3.0) * (x + y + z);
	double xr = r - x, yr = r - y, zr = r - z;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
	 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Z coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
	 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
	 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
	 */
double noise3_XYBeforeZ(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z){

	// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xy = x + y;
	double s2 = xy * -0.211324865405187;
	double zz = z * 0.577350269189626;
	double xr = x + s2 - zz, yr = y + s2 - zz;
	double zr = xy * 0.577350269189626 + zz;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
	 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Y coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
	 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
	 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
	 */
double noise3_XZBeforeY(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z){

	// Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xz = x + z;
	double s2 = xz * -0.211324865405187;
	double yy = y * 0.577350269189626;
	double xr = x + s2 - yy;
	double zr = z + s2 - yy;
	double yr = xz * 0.577350269189626 + yy;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}

/**
	 * 4D SuperSimplex noise base.
	 * Using ultra-simple 4x4x4x4 lookup partitioning.
	 * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
	 * but it does compete performance-wise with optimized OpenSimplex1.
	 */
double _noise4_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys, double zs, double ws){
	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs), ysb = fast_floor(ys), zsb = fast_floor(zs), wsb = fast_floor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// Unskewed offsets
	double ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011;
	double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;

	int index = ((fast_floor(xs * 4) & 3) << 0) | ((fast_floor(ys * 4) & 3) << 2) | ((fast_floor(zs * 4) & 3) << 4) | ((fast_floor(ws * 4) & 3) << 6);

	// Point contributions
	LatticePoint4D *c = (LatticePoint4D *)ose->LOOKUP_4D[index].data;
	for (int i = 0; i < ose->LOOKUP_4D[index].length; i++){
		double dx = xi + c->dx, dy = yi + c->dy, dz = zi + c->dz, dw = wi + c->dw;
		double attn = 0.8 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0){
			attn *= attn;

			int pxm = (xsb + c->xsv) & PMASK, pym = (ysb + c->ysv) & PMASK;
			int pzm = (zsb + c->zsv) & PMASK, pwm = (wsb + c->wsv) & PMASK;
			Grad4 grad = osg->permGrad4[osg->perm[osg->perm[osg->perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double extrapolation = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			value += attn * attn * extrapolation;
		}
		c++;
	}
	return value;
}

/**
	 * 4D SuperSimplex noise, classic lattice orientation.
	 */
double noise4_Classic(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w){

	// Get points for A4 lattice
	double s = 0.309016994374947 * (x + y + z + w);
	double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
	 * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
	 * Recommended for noise(x, y, sin(time), cos(time)) trick.
	 */
double noise4_XYBeforeZW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w){

	double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
	double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
	double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
	 * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
	 */
double noise4_XZBeforeYW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w){

	double s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718;
	double t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942;
	double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}

/**
	 * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
	 * and W for an extra degree of freedom.
	 * Recommended for time-varied animations which texture a 3D object (W=time)
	 */
double noise4_XYZBeforeW(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double x, double y, double z, double w){

	double xyz = x + y + z;
	double ww = w * 1.118033988749894;
	double s2 = xyz * -0.16666666666666666 + ww;
	double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}