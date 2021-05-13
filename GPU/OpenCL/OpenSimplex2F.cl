#define PSIZE 2048
#define PMASK 2047
#define PERIOD 64.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD





typedef struct {
    int xsv, ysv;
    double dx, dy;
} LatticePoint2D;

typedef struct {
    double dxr, dyr, dzr;
    int xrv, yrv, zrv;
	bool is_null;
} _LatticePoint3D;

typedef struct {
    _LatticePoint3D _this;
    _LatticePoint3D nextOnFailure;
    _LatticePoint3D nextOnSuccess;
} LatticePoint3D;

typedef struct {
    int xsv, ysv, zsv, wsv;
    double dx, dy, dz, dw;
    double xsi, ysi, zsi, wsi;
    double ssiDelta;
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

typedef struct {
    short perm[PSIZE];
    Grad2 permGrad2[PSIZE];
    Grad3 permGrad3[PSIZE];
    Grad4 permGrad4[PSIZE];
} OpenSimplexGradients;

typedef struct {
    Grad2 GRADIENTS_2D[PSIZE];
    Grad3 GRADIENTS_3D[PSIZE];
    Grad4 GRADIENTS_4D[PSIZE];
    LatticePoint2D LOOKUP_2D[4];
    LatticePoint3D LOOKUP_3D[8];
    LatticePoint4D VERTICES_4D[16];
} OpenSimplexEnv;



int get_index(){
	int x = get_global_id(0);
	int y = get_global_id(1);
	return y*get_global_size(0) + x;
}

double get_x(){
	int x = get_global_id(0);
	return (x + OFF_X) * FREQ;
}

double get_y(){
	int y = get_global_id(1);
	return (y + OFF_Y) * FREQ;
}

/*
	 * Utility
	 */

int _fastFloor(double x){
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
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
__kernel void noise2(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, const unsigned int size, __global double* output){

	int index = get_index();
    if (index < size){

		double x = get_x();
		double y = get_y();
		// Get points for A2* lattice
		double s = 0.366025403784439 * (x + y);
		double xs = x + s, ys = y + s;

		output[index] = _noise2_Base(ose, osg, xs, ys);
	}
}

/**
	 * 2D Simplex noise, with Y pointing down the main diagonal.
	 * Might be better for a 2D sandbox style game, where Y is vertical.
	 * Probably slightly less optimal for heightmaps or continent maps.
	 */
/*__kernel void noise2_XBeforeY(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y){

	// Skew transform and rotation baked into one.
	double xx = x * 0.7071067811865476;
	double yy = y * 1.224744871380249;

	return _noise2_Base(ose, osg, yy + xx, yy - xx);
}*/

/**
	 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
	 * Lookup table implementation inspired by DigitalShadow.
	 * It was actually faster to narrow down the points in the loop itself,
	 * than to build up the index with enough info to isolate 4 points.
	 */
/*double _noise3_BCC(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xr, double yr, double zr){

	// Get base and offsets inside cube of first lattice.
	int xrb = _fastFloor(xr), yrb = _fastFloor(yr), zrb = _fastFloor(zr);
	double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	_LatticePoint3D *c = &(ose->LOOKUP_3D[index]._this);
	while (!c->is_null){
		double dxr = xri + c->dxr, dyr = yri + c->dyr, dzr = zri + c->dzr;
		double attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0){
			c = &(ose->LOOKUP_3D[index].nextOnFailure);
		}
		else{
			int pxm = (xrb + c->xrv) & PMASK, pym = (yrb + c->yrv) & PMASK, pzm = (zrb + c->zrv) & PMASK;
			Grad3 grad = osg->permGrad3[osg->perm[osg->perm[pxm] ^ pym] ^ pzm];
			double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

			attn *= attn;
			value += attn * attn * extrapolation;
			c = &(ose->LOOKUP_3D[index].nextOnSuccess);
		}
	}
	return value;
}*/

/**
	 * 3D Re-oriented 4-point BCC noise, classic orientation.
	 * Proper substitute for 3D Simplex in light of Forbidden Formulae.
	 * Use noise3_XYBeforeZ or noise3_XZBeforeY instead, wherever appropriate.
	 */
/*__kernel void noise3_Classic(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z){

	// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
	// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
	// Orthonormal rotation. Not a skew transform.
	double r = (2.0 / 3.0) * (x + y + z);
	double xr = r - x, yr = r - y, zr = r - z;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}*/

/**
	 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Y).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Z coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
	 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
	 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
	 */
/*__kernel void noise3_XYBeforeZ(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z){

	// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
	// Orthonormal rotation. Not a skew transform.
	double xy = x + y;
	double s2 = xy * -0.211324865405187;
	double zz = z * 0.577350269189626;
	double xr = x + s2 - zz, yr = y + s2 - zz;
	double zr = xy * 0.577350269189626 + zz;

	// Evaluate both lattices to form a BCC lattice.
	return _noise3_BCC(ose, osg, xr, yr, zr);
}*/

/**
	 * 3D Re-oriented 4-point BCC noise, with better visual isotropy in (X, Z).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Y coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
	 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
	 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
	 */
/*__kernel void noise3_XZBeforeY(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z){

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
}*/

/**
	 * 4D OpenSimplex2F noise base.
	 * Current implementation not fully optimized by lookup tables.
	 * But still comes out slightly ahead of Gustavson's Simplex in tests.
	 */
/*double _noise4_Base(OpenSimplexEnv *ose, OpenSimplexGradients *osg, double xs, double ys, double zs, double ws){
	double value = 0;

	// Get base points and offsets
	int xsb = _fastFloor(xs), ysb = _fastFloor(ys), zsb = _fastFloor(zs), wsb = _fastFloor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
	double siSum = xsi + ysi + zsi + wsi;
	double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
	int inLowerHalf = (siSum < 2);
	if (inLowerHalf){
		xsi = 1 - xsi;
		ysi = 1 - ysi;
		zsi = 1 - zsi;
		wsi = 1 - wsi;
		siSum = 4 - siSum;
	}

	// Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
	double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
	double aabbScore = fabs(aabb), ababScore = fabs(abab), abbaScore = fabs(abba);

	// Find the closest point on the stretched tesseract as if it were the upper half
	int vertexIndex, via, vib;
	double asi, bsi;
	if (aabbScore > ababScore && aabbScore > abbaScore){
		if (aabb > 0){
			asi = zsi;
			bsi = wsi;
			vertexIndex = 0b0011;
			via = 0b0111;
			vib = 0b1011;
		}
		else{
			asi = xsi;
			bsi = ysi;
			vertexIndex = 0b1100;
			via = 0b1101;
			vib = 0b1110;
		}
	}
	else if (ababScore > abbaScore){
		if (abab > 0){
			asi = ysi;
			bsi = wsi;
			vertexIndex = 0b0101;
			via = 0b0111;
			vib = 0b1101;
		}
		else{
			asi = xsi;
			bsi = zsi;
			vertexIndex = 0b1010;
			via = 0b1011;
			vib = 0b1110;
		}
	}
	else{
		if (abba > 0){
			asi = ysi;
			bsi = zsi;
			vertexIndex = 0b1001;
			via = 0b1011;
			vib = 0b1101;
		}
		else{
			asi = xsi;
			bsi = wsi;
			vertexIndex = 0b0110;
			via = 0b0111;
			vib = 0b1110;
		}
	}
	if (bsi > asi){
		via = vib;
		double temp = bsi;
		bsi = asi;
		asi = temp;
	}
	if (siSum + asi > 3){
		vertexIndex = via;
		if (siSum + bsi > 4){
			vertexIndex = 0b1111;
		}
	}

	// Now flip back if we're actually in the lower half.
	if (inLowerHalf){
		xsi = 1 - xsi;
		ysi = 1 - ysi;
		zsi = 1 - zsi;
		wsi = 1 - wsi;
		vertexIndex ^= 0b1111;
	}

	// Five points to add, total, from five copies of the A4 lattice.
	for (int i = 0; i < 5; i++){

		// Update xsb/etc. and add the lattice point's contribution.
		LatticePoint4D *c = &(ose->VERTICES_4D[vertexIndex]);
		xsb += c->xsv;
		ysb += c->ysv;
		zsb += c->zsv;
		wsb += c->wsv;
		double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
		double dx = xi + c->dx, dy = yi + c->dy, dz = zi + c->dz, dw = wi + c->dw;
		double attn = 0.5 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0){
			int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
			Grad4 grad = osg->permGrad4[osg->perm[osg->perm[osg->perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double ramped = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			attn *= attn;
			value += attn * attn * ramped;
		}

		// Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
		if (i == 4)
			break;

		// Update the relative skewed coordinates to reference the vertex we just added.
		// Rather, reference its counterpart on the lattice copy that is shifted down by
		// the vector <-0.2, -0.2, -0.2, -0.2>
		xsi += c->xsi;
		ysi += c->ysi;
		zsi += c->zsi;
		wsi += c->wsi;
		ssi += c->ssiDelta;

		// Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
		double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
		vertexIndex = 0b0000;
		if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0){
			vertexIndex = 0b0001;
		}
		else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0){
			vertexIndex = 0b0010;
		}
		else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0){
			vertexIndex = 0b0100;
		}
		else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0){
			vertexIndex = 0b1000;
		}
	}

	return value;
}*/

/**
	 * 4D OpenSimplex2F noise, classic lattice orientation.
	 */
/*__kernel void noise4_Classic(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z, double w){

	// Get points for A4 lattice
	double s = -0.138196601125011 * (x + y + z + w);
	double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}*/

/**
	 * 4D OpenSimplex2F noise, with XY and ZW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
	 * Recommended for noise(x, y, sin(time), cos(time)) trick.
	 */
/*__kernel void noise4_XYBeforeZW(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z, double w){

	double s2 = (x + y) * -0.178275657951399372 + (z + w) * 0.215623393288842828;
	double t2 = (z + w) * -0.403949762580207112 + (x + y) * -0.375199083010075342;
	double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}*/

/**
	 * 4D OpenSimplex2F noise, with XZ and YW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
	 */
/*__kernel void noise4_XZBeforeYW(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z, double w){

	double s2 = (x + z) * -0.178275657951399372 + (y + w) * 0.215623393288842828;
	double t2 = (y + w) * -0.403949762580207112 + (x + z) * -0.375199083010075342;
	double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}*/

/**
	 * 4D OpenSimplex2F noise, with XYZ oriented like noise3_Classic,
	 * and W for an extra degree of freedom. W repeats eventually.
	 * Recommended for time-varied animations which texture a 3D object (W=time)
	 */
/*__kernel void noise4_XYZBeforeW(__global OpenSimplexEnv *ose, __global OpenSimplexGradients *osg, double x, double y, double z, double w){

	double xyz = x + y + z;
	double ww = w * 0.2236067977499788;
	double s2 = xyz * -0.16666666666666666 + ww;
	double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

	return _noise4_Base(ose, osg, xs, ys, zs, ws);
}*/
