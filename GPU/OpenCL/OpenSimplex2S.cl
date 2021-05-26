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
double _noise2_Base(short* perm, Grad2* permGrad2, LatticePoint2D* LOOKUP_2D, double xs, double ys){
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
__kernel void noise2(
	__global short* perm,
    __global Grad2* permGrad2,
    __global LatticePoint2D* LOOKUP_2D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		// Get points for A2* lattice
		double s = 0.366025403784439 * (x + y);
		double xs = x + s, ys = y + s;

		return _noise2_Base(ose, osg, xs, ys);
	}
}

/**
	 * 2D SuperSimplex noise, with Y pointing down the main diagonal.
	 * Might be better for a 2D sandbox style game, where Y is vertical.
	 * Probably slightly less optimal for heightmaps or continent maps.
	 */
__kernel void noise2_XBeforeY(
	__global short* perm,
    __global Grad2* permGrad2,
    __global LatticePoint2D* LOOKUP_2D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		// Skew transform and rotation baked into one.
		double xx = x * 0.7071067811865476;
		double yy = y * 1.224744871380249;

		return _noise2_Base(ose, osg, yy + xx, yy - xx);
	}
}

/**
	 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
	 * Lookup table implementation inspired by DigitalShadow.
	 * It was actually faster to narrow down the points in the loop itself,
	 * than to build up the index with enough info to isolate 8 points.
	 */
double _noise3_BCC(short* perm, Grad3* permGrad3, LatticePoint3D* LOOKUP_3D, double xr, double yr, double zr){

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
__kernel void noise3_Classic(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
		// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
		// Orthonormal rotation. Not a skew transform.
		double r = (2.0 / 3.0) * (x + y + z);
		double xr = r - x, yr = r - y, zr = r - z;

		// Evaluate both lattices to form a BCC lattice.
		return _noise3_BCC(ose, osg, xr, yr, zr);
	}
}

/**
	 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Y).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Z coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XYBeforeZ(x, z, Y) or use noise3_XZBeforeY.
	 * If Z is vertical in world coordinates, call noise3_XYBeforeZ(x, y, Z).
	 * For a time varied animation, call noise3_XYBeforeZ(x, y, T).
	 */
__kernel void noise3_XYBeforeZ(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

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
}

/**
	 * 3D Re-oriented 8-point BCC noise, with better visual isotropy in (X, Z).
	 * Recommended for 3D terrain and time-varied animations.
	 * The Y coordinate should always be the "different" coordinate in your use case.
	 * If Y is vertical in world coordinates, call noise3_XZBeforeY(x, Y, z).
	 * If Z is vertical in world coordinates, call noise3_XZBeforeY(x, Z, y) or use noise3_XYBeforeZ.
	 * For a time varied animation, call noise3_XZBeforeY(x, T, y) or use noise3_XYBeforeZ.
	 */
__kernel void noise3_XZBeforeY(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

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
}

/**
	 * 4D SuperSimplex noise base.
	 * Using ultra-simple 4x4x4x4 lookup partitioning.
	 * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
	 * but it does compete performance-wise with optimized OpenSimplex1.
	 */
double _noise4_Base(short* perm, Grad4* permGrad4, LatticePoint4D* VERTICES_4D, double xs, double ys, double zs, double ws){
	const unsigned int sizes[256] = {
		20, 15, 16, 17, 15, 16, 12, 15, 16, 12, 10, 14, 17, 15, 14,
		17, 15, 16, 12, 15, 16, 14, 14, 13, 12, 14, 11, 12, 15, 13,
		12, 14, 16, 12, 10, 14, 12, 14, 11, 12, 10, 11, 10, 13, 14,
		12, 13, 15, 17, 15, 14, 17, 15, 13, 12, 14, 14, 12, 13, 15,
		17, 14, 15, 17, 15, 16, 12, 15, 16, 14, 14, 13, 12, 14, 11,
		12, 15, 13, 12, 14, 16, 14, 14, 13, 14, 16, 16, 10, 14, 16,
		19, 11, 13, 10, 11, 10, 12, 14, 11, 12, 14, 16, 19, 11, 11,
		19, 13, 14, 12, 11, 14, 13, 15, 13, 12, 14, 13, 10, 11, 10,
		12, 11, 14, 13, 14, 10, 13, 13, 16, 12, 10, 14, 12, 14, 11,
		12, 10, 11, 10, 13, 14, 12, 13, 15, 12, 14, 11, 12, 14, 16,
		19, 11, 11, 19, 13, 14, 12, 11, 14, 13, 10, 11, 10, 13, 11,
		19, 13, 14, 10, 13, 16, 14, 13, 14, 14, 16, 14, 12, 13, 15,
		12, 11, 14, 13, 13, 14, 14, 16, 15, 13, 16, 15, 17, 15, 14,
		17, 15, 13, 12, 14, 14, 12, 13, 15, 17, 14, 15, 17, 15, 13,
		12, 14, 13, 10, 11, 10, 12, 11, 14, 13, 14, 10, 13, 13, 14,
		12, 13, 15, 12, 11, 14, 13, 13, 14, 14, 16, 15, 13, 16, 15,
		17, 14, 15, 17, 14, 10, 13, 13, 15, 13, 16, 15, 17, 13, 15, 20};

	const unsigned int offsets[256] = [0, 20, 35, 51, 68, 83, 99, 111, 126, 142, 154, 164, 178,
	195, 210, 224, 241, 256, 272, 284, 299, 315, 329, 343, 356, 368,
	382, 393, 405, 420, 433, 445, 459, 475, 487, 497, 511, 523, 537,
	548, 560, 570, 581, 591, 604, 618, 630, 643, 658, 675, 690, 704,
	721, 736, 749, 761, 775, 789, 801, 814, 829, 846, 860, 875, 892,
	907, 923, 935, 950, 966, 980, 994, 1007, 1019, 1033, 1044, 1056,
	1071, 1084, 1096, 1110, 1126, 1140, 1154, 1167, 1181, 1197, 1213,
	1223, 1237, 1253, 1272, 1283, 1296, 1306, 1317, 1327, 1339, 1353,
	1364, 1376, 1390, 1406, 1425, 1436, 1447, 1466, 1479, 1493, 1505,
	1516, 1530, 1543, 1558, 1571, 1583, 1597, 1610, 1620, 1631, 1641,
	1653, 1664, 1678, 1691, 1705, 1715, 1728, 1741, 1757, 1769, 1779,
	1793, 1805, 1819, 1830, 1842, 1852, 1863, 1873, 1886, 1900, 1912,
	1925, 1940, 1952, 1966, 1977, 1989, 2003, 2019, 2038, 2049, 2060,
	2079, 2092, 2106, 2118, 2129, 2143, 2156, 2166, 2177, 2187, 2200,
	2211, 2230, 2243, 2257, 2267, 2280, 2296, 2310, 2323, 2337, 2351,
	2367, 2381, 2393, 2406, 2421, 2433, 2444, 2458, 2471, 2484, 2498,
	2512, 2528, 2543, 2556, 2572, 2587, 2604, 2619, 2633, 2650, 2665,
	2678, 2690, 2704, 2718, 2730, 2743, 2758, 2775, 2789, 2804, 2821,
	2836, 2849, 2861, 2875, 2888, 2898, 2909, 2919, 2931, 2942, 2956,
	2969, 2983, 2993, 3006, 3019, 3033, 3045, 3058, 3073, 3085, 3096,
	3110, 3123, 3136, 3150, 3164, 3180, 3195, 3208, 3224, 3239, 3256,
	3270, 3285, 3302, 3316, 3326, 3339, 3352, 3367, 3380, 3396, 3411,
	3428, 3441, 3456];

	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs), ysb = fast_floor(ys), zsb = fast_floor(zs), wsb = fast_floor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// Unskewed offsets
	double ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011;
	double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;

	int index = ((fast_floor(xs * 4) & 3) << 0) | ((fast_floor(ys * 4) & 3) << 2) | ((fast_floor(zs * 4) & 3) << 4) | ((fast_floor(ws * 4) & 3) << 6);

	// Point contributions
	LatticePoint4D* c = &(ose->LOOKUP_4D[offsets[index]]);
	for (int i = 0; i < sizes[index]; i++){
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
__kernel void noise4_Classic(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		// Get points for A4 lattice
		double s = 0.309016994374947 * (x + y + z + w);
		double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

		return _noise4_Base(ose, osg, xs, ys, zs, ws);
	}
}

/**
	 * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
	 * Recommended for noise(x, y, sin(time), cos(time)) trick.
	 */
__kernel void noise4_XYBeforeZW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
		double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
		double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

		return _noise4_Base(ose, osg, xs, ys, zs, ws);
	}
}

/**
	 * 4D SuperSimplex noise, with XZ and YW forming orthogonal triangular-based planes.
	 * Recommended for 3D terrain, where X and Z (or Y and W) are horizontal.
	 */
__kernel void noise4_XZBeforeYW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		double s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718;
		double t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942;
		double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

		return _noise4_Base(ose, osg, xs, ys, zs, ws);
	}
}

/**
	 * 4D SuperSimplex noise, with XYZ oriented like noise3_Classic,
	 * and W for an extra degree of freedom.
	 * Recommended for time-varied animations which texture a 3D object (W=time)
	 */
__kernel void noise4_XYZBeforeW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	const unsigned int width,
	const unsigned int height,
	__global double* output){

	int index = get_index(width, height);
    if (index >= 0){

		double x = get_x();
		double y = get_y();

		double xyz = x + y + z;
		double ww = w * 1.118033988749894;
		double s2 = xyz * -0.16666666666666666 + ww;
		double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

		return _noise4_Base(ose, osg, xs, ys, zs, ws);
	}
}